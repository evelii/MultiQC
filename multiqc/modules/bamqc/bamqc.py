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

        # FInd and load any BamQC json files
        self.bamqc_data = dict()
        for f in self.find_log_files('bamqc'):
            # remove all file extensions
            file_name = f['s_name'].replace(".annotated", "")
            self.parse_bamqc_log(f['f'], file_name, f)

        # check reads for each sample and mark any outlier as 'fail'; otherwise, mark it as 'pass'    
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

        # Basic Stats Table
#       self.bamqc_general_table()

        # Add the statuses to the intro for multiqc_bamqc.js to pick up
        statuses = dict()
        for s_name in self.bamqc_data:
            for section, status in self.bamqc_data[s_name]['statuses'].items():
                    section = section.lower().replace(' ', '_')
                    try:
                            statuses[section][s_name] = status
                    except KeyError:
                            statuses[section] = {s_name: status}
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
                        'mapped reads', 'aligned bases', 'soft clip bases', 'insert mean']

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

                # reset the list
                reads = [] 

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
                'color': '#3ea540'  # passed: bars are plotted in green
        }
        cats['reads on target fail'] = {
                'name': 'Reads on Target (Outliers)',
                'color': '#bf3939'  # failed: bars are plotted in red
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
                'color': '#3ea540'  # green
        }
        cats['reads per start point fail'] = {
                'name': 'Reads per Start Point (Outliers)',
                'color': '#bf3939'  # red
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
                'color': '#3ea540'   # green
        }
        cats['total reads fail'] = {
                'name': 'Total Reads (Outliers)',
                'color': '#bf3939'   # red
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
                'color': '#3ea540'   # green
        }
        cats['paired reads fail'] = {
                'name': 'Paired Reads (Outliers)',
                'color': '#bf3939'    # red
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
                'color': '#3ea540'    # green
        }
        cats['mapped reads fail'] = {
                'name': 'Mapped Reads (Outliers)',
                'color': '#bf3939'    # red
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
                'color': '#3ea540'    # green
        }
        cats['aligned bases fail'] = {
                'name': 'Aligned Bases (Outliers)',
                'color': '#bf3939'    # red
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
                'color': '#3ea540'    # green
        }
        cats['soft clip bases fail'] = {
                'name': 'Soft Clip Bases (Outliers)',
                'color': '#bf3939'    # red
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
                'color': '#3ea540'   # green
        }
        cats['insert mean fail'] = {
                'name': 'Insert Mean (Outliers)',
                'color': '#bf3939'   # red
        }

        self.add_section (
                name = 'Insert Mean',
                anchor = 'bamqc_insert_mean',
                description = 'Insert mean data for each sample.',
                plot = bargraph.plot(bdata, cats, pconfig)
        ) 
