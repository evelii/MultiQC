#!/usr/bin/env python

""" MultiQC module to parse output from BamQC """

from __future__ import print_function
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule
import logging
import os
import io
import json

# Initialise the logger
log = logging.getLogger(__name__)

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

        # Filter to strip out ignored sample names
        self.bamqc_data = self.ignore_samples(self.bamqc_data)

        if len(self.bamqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bamqc_data)))

        # Write parsed report data to a file 
        self.write_data_file(self.bamqc_data, 'multiqc_bamqc')

        # Basic Stats Table
#       self.bamqc_general_table()

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
                bdata[sample]['reads on target'] = self.bamqc_data[sample]['reads on target']

        self.add_section (
                name = 'Reads on Target',
                anchor = 'bamqc_reads_on_target',
                description = 'Reads on target data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['reads per start point'] = self.bamqc_data[sample]['reads per start point']

        self.add_section (
                name = 'Reads per Start Point',
                anchor = 'bamqc_reads_per_start_point',
                description = 'Reads per Start Point data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['total reads'] = self.bamqc_data[sample]['total reads']

        self.add_section (
                name = 'Total Reads',
                anchor = 'bamqc_total_reads',
                description = 'Total Reads data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['paired reads'] = self.bamqc_data[sample]['paired reads']

        self.add_section (
                name = 'Paired Reads',
                anchor = 'bamqc_paired_reads',
                description = 'Paired reads data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['mapped reads'] = self.bamqc_data[sample]['mapped reads']

        self.add_section (
                name = 'Mapped Reads',
                anchor = 'bamqc_mapped_reads',
                description = 'Mapped Reads data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['aligned bases'] = self.bamqc_data[sample]['aligned bases']

        self.add_section (
                name = 'Aligned Bases',
                anchor = 'bamqc_aligned_bases',
                description = 'Aligned Bases data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['soft clip bases'] = self.bamqc_data[sample]['soft clip bases']

        self.add_section (
                name = 'Soft Clip Bases',
                anchor = 'bamqc_soft_clip_bases',
                description = 'Soft clip bases data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
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
                bdata[sample]['insert mean'] = self.bamqc_data[sample]['insert mean']

        self.add_section (
                name = 'Insert Mean',
                anchor = 'bamqc_insert_mean',
                description = 'Insert mean data for each sample.',
                plot = bargraph.plot(bdata, None, pconfig)
        ) 
