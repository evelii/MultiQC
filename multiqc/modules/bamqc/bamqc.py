#!/usr/bin/env python

""" MultiQC module to parse output from BamQC """
from __future__ import print_function
import sys
import logging
# Initialise the logger
log = logging.getLogger(__name__)

from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph
from multiqc.plots import table
from multiqc.utils import report
from multiqc.modules.base_module import BaseMultiqcModule
import os
import io
import json
import statistics

class MultiqcModule(BaseMultiqcModule):
    """ BamQC module """ 

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='BamQC', anchor='bamqc', 
        href="https://github.com/oicr-gsi/BamQC", 
        info="is an application for analysing BAM files containing mapped data and generating a QC report.")

        # Find and load any BamQC json files
        self.bamqc_data = dict()
        # Store short sample Names, Run Names and Indices
        self.sample_name_info = dict() 

        for f in self.find_log_files('bamqc'):
            # remove all file extensions
            file_name = f['s_name'].replace(".annotated", "")
            self.parse_bamqc_log(f['f'], file_name, f)

        # check reads for each sample and mark any outlier as 'fail'; otherwise, mark it as 'pass'    
        # bamqc_stats_warning is used to save whether a set of data has (mean - 2 standard deviations) < 0; if so, 
        # add a warning to the table near the progress bar
        self.bamqc_stats_warning = dict()
        self.find_outliers()

        # Filter to strip out ignored sample names
        self.bamqc_data = self.ignore_samples(self.bamqc_data)

        if len(self.bamqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bamqc_data)))

        # Write parsed report data to a file 
        self.write_data_file(self.bamqc_data, 'multiqc_bamqc')

        # Add to self.css and self.js to be included in template
        self.css = { 'assets/css/multiqc_bamqc.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_bamqc.css') }
        self.js = { 'assets/js/multiqc_bamqc.js' : os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_bamqc.js') }

        # BamQC Stats Table
        self.bamqc_general_stats()

        # Add the statuses to the intro for multiqc_bamqc.js to pick up
        statuses = dict()
        for s_name in self.bamqc_data:
            for section, status in self.bamqc_data[s_name]['statuses'].items():
                    section = section.lower().replace(' ', '_')
                    try:
                            statuses[section][s_name] = status
                    except KeyError:
                            statuses[section] = {s_name: status}
        for section in statuses.keys():
            statuses[section]['warning'] = self.bamqc_stats_warning[section]

        self.intro += '<script type="text/javascript">bamqc_passfails = {};</script>'.format(json.dumps(statuses))

        # Add each Sample Read Plot section in order
        self.reads_on_target_plot()
        self.reads_per_start_point_plot()
        self.total_reads_plot()
        self.paired_reads_plot()
        self.mapped_reads_plot()
        self.aligned_bases_plot()
        self.soft_clip_bases_plot()
        self.insert_mean_plot()

    def parse_bamqc_log(self, file_content, s_name, f):
        """ Takes contents from a .bam.BamQC.json file and parses out required
        statistics and data. Returns a dict with keys 'sample_name' and values
        in the form of a dict (e.g. 'reads on target: 34242'). """

        # s_name is the Sample name (from cleaned filename)
        if s_name in self.bamqc_data.keys():
           log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.add_data_source(f, s_name)
        # the value of a sample is a dictionary
        self.bamqc_data[s_name] = dict() 

        # properties that need to be parsed and analyzed
        property_list = ['reads on target', 'reads per start point', 'total reads', 'paired reads', \
                        'mapped reads', 'aligned bases', 'soft clip bases', 'insert mean', 'target size']

        self.sample_name_info[s_name] = dict()

        for s in file_content.split(','):
                property = s.split(':')[0]

                # remove '{' and '"' from property names
                property = property.replace("{", "")
                property = property.replace("\"", "")

                if property in property_list: 

                        if property in self.bamqc_data[s_name].keys():
                               log.debug("Duplicate property! Overwriting: {}".format(property))

                        # remove quotes from data values
                        value = s.split(':')[1].replace("\"", "")

                        # convert strings to numbers
                        if "." in value:
                               self.bamqc_data[s_name][str(property)] = float(value)
                        else:
                               self.bamqc_data[s_name][str(property)] = int(value)

                elif property == 'library':
                       
                        library = s.split(':')[1].replace("\"", "")
                        self.sample_name_info[s_name]['library'] = library

                elif property == 'run name':
 
                        run_name = s.split(':')[1].replace("\"", "")
                        self.sample_name_info[s_name]['run name'] = run_name

                elif property == 'barcode':

                        index = s.split(':')[1].replace("\"", "")
                        self.sample_name_info[s_name]['index'] = index

                elif property == 'group id':
 
                        group_id = s.split(':')[1].replace("\"", "")
                        self.sample_name_info[s_name]['group id'] = group_id

                elif property == 'lane':

                        lane_number = s.split(':')[1]
                        self.sample_name_info[s_name]['lane'] = lane_number

        # check if any required data is missing
        for p in property_list:
               if not (p in self.bamqc_data[s_name].keys()):
                       log.warning("Data for {} is missing!".format(p))

    # an outlier is a value which is more than 2 standard deviations away from the mean in a dataset
    def find_outliers(self):
        reads = []
        mean = 0 
        stdev = 0
        property_list = ['reads on target', 'reads per start point', 'total reads', 'paired reads', \
                        'mapped reads', 'aligned bases', 'soft clip bases', 'insert mean']
        for property in property_list:

                for s_name in self.bamqc_data.keys():
                        reads.append(self.bamqc_data[s_name][property])
                        # keeping track of if any read is an outlier (more than two standard deviations away from the mean)
                        if 'statuses' not in self.bamqc_data[s_name]:
                                    self.bamqc_data[s_name]['statuses'] = dict()

                try:
                        mean = statistics.mean(reads)
                        stdev = statistics.stdev(reads)
                except Exception as e:
                        log.warning("Not enough data is available to calculate mean and stdev!")
                        log.debug("Reads are '{}'".format(reads))

                for s_name in self.bamqc_data.keys():
                        s_reads = self.bamqc_data[s_name][property]
                        if s_reads > (mean + 2 * stdev) or s_reads < (mean - 2 * stdev):
                                   self.bamqc_data[s_name]['statuses'][property] = 'fail'
                        else:
                                   self.bamqc_data[s_name]['statuses'][property] = 'pass'
                        
                # Check if (mean - 2 * stdev) < 0; if so, set the value in bamqc_stats_warning to be true
                n_property = property.lower().replace(' ', '_')
                self.bamqc_stats_warning[n_property] = (mean - 2 * stdev) < 0

                # reset the list
                reads = [] 

    def bamqc_general_stats(self):
        """ Add single-number stats to the bamqc general statistics
        table at the top of the report but below the general stats table"""

        # Prepare the data
        
        # Calculate map %: 
        # formula = (mapped reads / total reads) * 100
        # if total reads is 0, treat as 1
        for s_name in self.bamqc_data.keys():
                mapped_reads = self.bamqc_data[s_name]['mapped reads']
                total_reads = self.bamqc_data[s_name]['total reads']
                if total_reads == 0:
                        total_reads = 1
                map_percent = (mapped_reads / total_reads) * 100
                self.bamqc_data[s_name]['map percent'] = map_percent

        # Calculate on target %:
        # formula = (reads on target / mapped reads) * 100
        # if mapped reads is 0, treat as 1
        for s_name in self.bamqc_data.keys():
                mapped_reads = self.bamqc_data[s_name]['mapped reads'] or 1
                reads_on_target = self.bamqc_data[s_name]['reads on target']
                on_target_percent = (reads_on_target / mapped_reads) * 100
                self.bamqc_data[s_name]['on target percent'] = on_target_percent

        # Calculate the estimated total coverage
        # formula = estimated yield / target size
        # if target size is 0, treat as 1        
        # formula for the estimated yield = (total aligned bases * on target percentage) / reads per start point
        # if reads per start point is 0, treat as 1
        for values in self.bamqc_data.values():
               target_size = values['target size'] or 1
               aligned_bases = values['aligned bases']
               on_target_percent = values['on target percent'] / 100
               reads_per_start_point = values['reads per start point'] or 1
               est_yield = (aligned_bases * on_target_percent) / reads_per_start_point
               coverage = est_yield / target_size
               values['coverage'] = coverage

        # Manipulate/Shorten Sample names 
        # Keep the date and the lane number in Sample Name
        # Move Run Name and Index into separate columns
        stats_data = dict(self.bamqc_data) # Make a copy of self.bamqc_data to contain new sample names and use it to construct the table

        for s_name in self.bamqc_data.keys():
                date = self.sample_name_info[s_name]['run name'].split('_')[0]
                lane_number = "L" + str(self.sample_name_info[s_name]['lane']).zfill(3)  # prepend zeros on the left to fill width of 3
                if 'group id' in self.sample_name_info[s_name].keys():
                        sample_name = self.sample_name_info[s_name]['library'] + "_" + self.sample_name_info[s_name]['group id'] + "_" + date + "_" + lane_number
                else:
                        sample_name = self.sample_name_info[s_name]['library'] + "_" + date + "_" + lane_number
                stats_data[sample_name] = stats_data[s_name]
                del stats_data[s_name]
                stats_data[sample_name]['run name'] = self.sample_name_info[s_name]['run name']
                stats_data[sample_name]['index'] = self.sample_name_info[s_name]['index']

        table_config = {
                'namespace': 'BamQC',
                'table_title': 'Aligned Statistics'
        }
 
        headers = OrderedDict()
        headers['run name'] = {
                'title': 'Run Name',
                'description': 'Run Name',
                'hidden': True
        }
        headers['index'] = {
                'title': 'Index',
                'description': 'Barcode/Index',
                'hidden': True
        }
        headers['map percent'] = {
                'title': 'Map %',
                'description': 'Map Percent',
                'scale': 'Blues'
        }
        headers['on target percent'] = {
                'title': '% On Target',
                'description': 'Percent on Target',
                'scale': 'Blues'
        }
        headers['coverage'] = {
                'title': 'Coverage',
                'description': 'Estimated Coverage',
                'scale': 'Blues'
        }
        headers['reads on target'] = {
                'title': 'On Target',
                'description': 'Reads on Target',
                'scale': 'Blues',
                'format': '{:,.0f}'  # no decimal places
        }
        headers['reads per start point'] = {
                'title': 'Reads/SP',
                'description': 'Reads per Start Point',
                'scale': 'Blues',
                'format': '{:.2f}'   # show two decimal places
        }
        headers['total reads'] = {
                'title': 'Total Rs',
                'description': 'Total Reads',
                'scale': 'Blues',
                'format': '{:,.0f}'
        }
        headers['paired reads'] = {
                'title': 'Paired Rs',
                'description': 'Paired Reads',
                'scale': 'Blues',
                'format': '{:,.0f}'
        }
        headers['mapped reads'] = {
                'title': 'Mapped Rs',
                'description': 'Mapped Reads',
                'scale': 'Blues',
                'format': '{:,.0f}'
        }
        headers['aligned bases'] = {
                'title': 'Aligned Rs',
                'description': 'Aligned Reads',
                'scale': 'Blues',
                'format': '{:,.0f}'
        }
        headers['soft clip bases'] = {
                'title': 'Soft Clip',
                'description': 'Soft Clip Bases',
                'scale': 'Blues',
                'format': '{:,.0f}'
        } 
        headers['insert mean'] = {
                'title': 'Insert Mean',
                'description': 'Insert Mean',
                'scale': 'Blues',
                # one decimal place will be used
        }
        table_html = table.plot(stats_data, headers, table_config)
        report.bamqc_general_stats_html = table_html  # add the table to a report

    def reads_on_target_plot(self):
#       Bar plots showing sample names and values for reads on target
        pconfig = {
                'id': 'bamqc_reads_on_target_plot',
                'title': 'BamQC: Reads on Target',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "reads on target" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['reads on target'] == 'fail':
                        bdata[sample]['reads on target fail'] = self.bamqc_data[sample]['reads on target']
                else:
                        bdata[sample]['reads on target'] = self.bamqc_data[sample]['reads on target']

        # specify colors to distinguish passed values and outliers
        cats = OrderedDict()
        cats['reads on target'] = {
                'name': 'Reads on Target',
                'color': '#55b4d1'  # passed: bars are plotted in blue
        }
        cats['reads on target fail'] = {
                'name': 'Reads on Target (Outliers)',
                'color': '#ddac39'  # failed: bars are plotted in yellow
        }

        self.add_section (
                name = 'Reads on Target',
                anchor = 'bamqc_reads_on_target',
                description = 'Reads on target data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def reads_per_start_point_plot(self):
#       Bar plots showing sample names and values for reads per start point
        pconfig = {
                'id': 'bamqc_reads_per_start_point_plot',
                'title': 'BamQC: Reads per Start Point',
                'cpswitch': False,
                'tt_percentages': False,
                'yDecimals': True,
                'tt_decimals': 2 
        }

        bdata = dict()
        # pull out information about "reads per start point" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['reads per start point'] == 'fail':
                        bdata[sample]['reads per start point fail'] = self.bamqc_data[sample]['reads per start point']
                else:
                        bdata[sample]['reads per start point'] = self.bamqc_data[sample]['reads per start point']

        cats = OrderedDict()
        cats['reads per start point'] = {
                'name': 'Reads per Start Point',
                'color': '#55b4d1'  # blue
        }
        cats['reads per start point fail'] = {
                'name': 'Reads per Start Point (Outliers)',
                'color': '#ddac39'  # yellow
        }

        self.add_section (
                name = 'Reads per Start Point',
                anchor = 'bamqc_reads_per_start_point',
                description = 'Reads per Start Point data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def total_reads_plot(self):
#       Bar plots showing sample names and values for total reads 
        pconfig = {
                'id': 'bamqc_total_reads_plot',
                'title': 'BamQC: Total Reads',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "total reads" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['total reads'] == 'fail':
                        bdata[sample]['total reads fail'] = self.bamqc_data[sample]['total reads']
                else:
                        bdata[sample]['total reads'] = self.bamqc_data[sample]['total reads']

        cats = OrderedDict()
        cats['total reads'] = {
                'name': 'Total Reads',
                'color': '#55b4d1'   # blue
        }
        cats['total reads fail'] = {
                'name': 'Total Reads (Outliers)',
                'color': '#ddac39'   # yellow
        }

        self.add_section (
                name = 'Total Reads',
                anchor = 'bamqc_total_reads',
                description = 'Total Reads data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def paired_reads_plot(self):
#       Bar plots showing sample names and paired reads 
        pconfig = {
                'id': 'bamqc_paired_reads_plot',
                'title': 'BamQC: Paired Reads',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "paired reads" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['paired reads'] == 'fail':
                        bdata[sample]['paired reads fail'] = self.bamqc_data[sample]['paired reads']
                else:
                        bdata[sample]['paired reads'] = self.bamqc_data[sample]['paired reads']

        cats = OrderedDict()
        cats['paired reads'] = {
                'name': 'Paired Reads',
                'color': '#55b4d1'   # blue
        }
        cats['paired reads fail'] = {
                'name': 'Paired Reads (Outliers)',
                'color': '#ddac39'    # yellow
        }

        self.add_section (
                name = 'Paired Reads',
                anchor = 'bamqc_paired_reads',
                description = 'Paired reads data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def mapped_reads_plot(self):
#       Bar plots showing sample names and mapped reads
        pconfig = {
                'id': 'bamqc_mapped_reads_plot',
                'title': 'BamQC: Mapped Reads',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "mapped reads" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['mapped reads'] == 'fail':
                        bdata[sample]['mapped reads fail'] = self.bamqc_data[sample]['mapped reads']
                else:
                        bdata[sample]['mapped reads'] = self.bamqc_data[sample]['mapped reads']

        cats = OrderedDict() 
        cats['mapped reads'] = {
                'name': 'Mapped Reads',
                'color': '#55b4d1'    # blue
        }
        cats['mapped reads fail'] = {
                'name': 'Mapped Reads (Outliers)',
                'color': '#ddac39'    # yellow
        }

        self.add_section (
                name = 'Mapped Reads',
                anchor = 'bamqc_mapped_reads',
                description = 'Mapped Reads data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def aligned_bases_plot(self):
#       Bar plots showing sample names and aligned bases 
        pconfig = {
                'id': 'bamqc_aligned_bases_plot',
                'title': 'BamQC: Aligned Bases',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "aligned bases" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['aligned bases'] == 'fail':
                        bdata[sample]['aligned bases fail'] = self.bamqc_data[sample]['aligned bases']
                else:
                        bdata[sample]['aligned bases'] = self.bamqc_data[sample]['aligned bases']

        cats = OrderedDict()
        cats['aligned bases'] = {
                'name': 'Aligned Bases',
                'color': '#55b4d1'    # blue
        }
        cats['aligned bases fail'] = {
                'name': 'Aligned Bases (Outliers)',
                'color': '#ddac39'    # yellow
        }

        self.add_section (
                name = 'Aligned Bases',
                anchor = 'bamqc_aligned_bases',
                description = 'Aligned Bases data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def soft_clip_bases_plot(self):
#       Bar plots showing sample names and soft clip bases 
        pconfig = {
                'id': 'bamqc_soft_clip_bases_plot',
                'title': 'BamQC: Soft Clip Bases',
                'cpswitch': False,
                'tt_percentages': False
        }

        bdata = dict()
        # pull out information about "soft clip bases" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['soft clip bases'] == 'fail':
                        bdata[sample]['soft clip bases fail'] = self.bamqc_data[sample]['soft clip bases']
                else:
                        bdata[sample]['soft clip bases'] = self.bamqc_data[sample]['soft clip bases']

        cats = OrderedDict()
        cats['soft clip bases'] = {
                'name': 'Soft Clip Bases',
                'color': '#55b4d1'    # blue
        }
        cats['soft clip bases fail'] = {
                'name': 'Soft Clip Bases (Outliers)',
                'color': '#ddac39'    # yellow
        }

        self.add_section (
                name = 'Soft Clip Bases',
                anchor = 'bamqc_soft_clip_bases',
                description = 'Soft clip bases data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        )

    def insert_mean_plot(self):
#       Bar plots showing sample names and insert mean 
        pconfig = {
                'id': 'bamqc_insert_mean_plot',
                'title': 'BamQC: Insert Mean',
                'cpswitch': False,
                'tt_percentages': False,
                'yDecimals': True,
                'tt_decimals': 1
        }

        bdata = dict()
        # pull out information about "insert mean" from bamqc_data dictionary
        for sample in self.bamqc_data.keys():
                bdata[sample] = dict()
                if self.bamqc_data[sample]['statuses']['insert mean'] == 'fail':
                        bdata[sample]['insert mean fail'] = self.bamqc_data[sample]['insert mean']
                else:
                        bdata[sample]['insert mean'] = self.bamqc_data[sample]['insert mean']

        cats = OrderedDict()
        cats['insert mean'] = {
                'name': 'Insert Mean',
                'color': '#55b4d1'   # blue
        }
        cats['insert mean fail'] = {
                'name': 'Insert Mean (Outliers)',
                'color': '#ddac39'   # yellow
        }

        self.add_section (
                name = 'Insert Mean',
                anchor = 'bamqc_insert_mean',
                description = 'Insert mean data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        ) 
