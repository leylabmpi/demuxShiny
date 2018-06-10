#!/usr/bin/env python
from __future__ import print_function
import sys,os
import argparse
from pprint import pprint
from datetime import datetime
from collections import OrderedDict, Counter
from schema import Schema, And, Or, Use, Optional, Regex, SchemaError

#-- argparse--#
desc = 'Validate demultiplexing pipeline samples sheet format'
epi = """DESCRIPTION:
  See https://confluence.eb.local:8443/display/M005SEQ/Demultiplexing+pipeline+for+Illumina+sequencers
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('samples_sheet', metavar='samples_sheet', type=str,
                    help='Samples sheet file path')
parser.add_argument('-s', '--seq-tech', type=str, default='HiSeq',
                    choices=['HiSeq', 'MiSeq'],
                    help='Sequencing technology (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')


class Parser(object):
    """Class for parsing Samples Sheet sections
    """
    def __init__(self):
        self.header = None
    
    def parse_params(self, x, line):
        """Parsing parameter-formatted sections
        """
        x = x.split(',')
        if len(x) < 2:
            x.append('')
        x = {x[0] : x[1]}
        return x
        
    def parse_table(self, x, line):
        """Parsing table-formatted sections.
        WARNING: currently can only handle 1 table
        """
        x = x.split(',')
        if self.header is None:
            self.header = x
        else:
            z = {}
            for i,y in enumerate(x):
                try:
                    k = self.header[i]
                except IndexError:
                    msg = 'ERROR in Line {}: row values don\'t match header'
                    print(msg.format(line))
                    sys.exit(1)
                z[k] = y
            return z
        return None

    def warn_extra(self, x, line):
        """Warning that there is content in the '[Extra]' key
        """
        x = x.split(',')
        msg = 'WARNING: extra content detected at line {}: {}'.format(line,x)
        print(msg)
    
def read_samples_sheet(path):
    """Reading in the Sample Sheet.
    path : path to Sample Sheet file
    """
    P = Parser()
    sections = {'[Header]' : P.parse_params,
                '[Reads]' : P.parse_params,
                '[Settings]' : P.parse_params,
                '[Data]' : P.parse_table,
                '[Extra]' : P.warn_extra}
    
    sheet = OrderedDict()
    with open(path) as inF:
        section = '[Extra]'
        for i,line in enumerate(inF):
            line = line.strip().rstrip(',')
            if line == '':
                continue
            elif line in sections:
                section = line                
            else:
                if line.startswith('['):
                    msg = 'WARNING in Line {}: line starts with "["'
                    msg += ', but {} is not a valid section ID'
                    print(msg.format(i, line))
                try:
                    line = sections[section](line, i)
                except KeyError:
                    msg = 'ERROR in Line {}: "{}" not a valid section'
                    print(msg.format(i))
                    sys.exit(1)
                if line is None:
                    continue
                try:
                    sheet[section].append(line)
                except KeyError:
                    sheet[section] = [line]

    return(sheet)

def header_schema():
    """Validation schema for header section
    """
    s = [{Optional('IEMFileVersion') : '4'},
         {Optional('Experiment Name') : str},
         {Optional('Date') : Regex(r'\d{1,2}\/\d{1,2}\/20\d{2}')},
         {'Workflow': Or('GenerateFASTQ', 'GenerateChromiumFASTQ')},
         {Optional('Application') : str},
         {Optional('Assay') : str},
         {Optional('Description') : str},
         {Optional('Chemistry') : str}]    
    return(Schema(s))

def reads_schema():
    """Validation schema for reader section
    """
    s = [{And(Use(int), lambda n: 0 < n <= 300) : ''}]
    msg = 'Read length must be > 0 and <= 300'    
    return(Schema(s, error=msg))

def settings_schema():
    """Validation schema for settings section
    """
    rgx = Regex(r'^[A-Z]+$')
    s = [{Optional(Or('Adapter', 'TrimAdapter ')) : rgx},
         {Optional(Or('AdapterRead2', 'TrimAdapterRead2')) : rgx},
         {Optional('MaskAdapter') : rgx},
         {Optional('MaskAdapterRead2') : rgx},
         {Optional('ReverseComplement') : Or('0','1')},
         {Optional('FindAdaptersWithIndels') : Or('0','1')}]         
    return(Schema(s))

def data_schema(seq_tech):
    """Validation schema for data section
    seq_tech : HiSeq or MiSeq
    """
    rgx = Regex(r'^[A-Za-z0-9_-]{1,100}$')
    idx_rgx = Regex(r'^[A-Z]+$')
    well_rgx = Regex(r'^[A-Z]0*[0-9]+$')
    s = {'Sample_ID' : rgx,
         'index' : idx_rgx,
         'Sample_Project' : rgx,
         'Description' : rgx,
         Optional('Sample_Name') : Or('', rgx),
         Optional('Sample_Plate') : Or('', rgx),
         Optional('Sample_Well') : Or('', well_rgx),
         Optional('I7_Index_ID') : Or('', rgx),
         Optional('I5_Index_ID') : Or('', rgx),
         Optional('index2') : Or('', idx_rgx)}
    if seq_tech.lower() == 'hiseq':
        s['Lane'] = Regex(r'^[1-8]$')
    return(Schema([s]))

def unique_values(rows, idx):
    """Checking for unique values in table.
    rows : [{ }, { }, ...]
    idx : dict key (column name)
    """
    # counting occurances of the values (value for "idx" key)
    cnt = {}
    for i,x in enumerate(rows):
        try:
            value = x[idx]
        except KeyError:
            msg = 'ERROR in table row {}: "{}" key not found'
            print(msg.format(i, idx))
        try:
            cnt[value] += 1
        except KeyError:
            cnt[value] = 1

    # error if duplicates found
    errors = []
    for k,v in cnt.items():
        if v > 1:
            msg = 'ERROR: "{}" is found {} times in table'
            errors.append(msg.format(k,v))
    return errors
            
def validate_by_line(section, f):
    errors = []
    for x in section:
        try:
            f.validate([x])
        except SchemaError as e:
            errors.append(e)
    return errors
            
def validate_samples_sheet(sheet, seq_tech):
    """Validating each section of the Samples Sheet
    sheet : parsed sheet object
    seq_tech : str; HiSeq, MiSeq, etc
    """
    # validation schemas for each section
    schemas = {'[Header]' : header_schema(),
               '[Reads]' : reads_schema(),
               '[Settings]' : settings_schema(),
               '[Data]' : data_schema(seq_tech)}

    # validating
    errors = {}
    for x,f in schemas.items():
        # checking that sections exist ([Settings] is optional)
        try:
            sheet[x]
        except KeyError:
            if x == '[Settings]':
                continue
            else:
                msg = 'ERROR: Cannot find section "{}"'
                print(msg.format(x))
                sys.exit(1)
        # validating
        errors[x] = validate_by_line(sheet[x], f)

    # checking for unique IDs in [Data]
    e = unique_values(sheet['[Data]'], 'Sample_ID')
    if len(e) > 0:
        errors['[Data]'].append(e)

    return errors

def validate_path(path):
    """Validate path
    path : str; file path
    """
    if not os.path.exists(path):
        msg = 'ERROR: cannot find file: {}'
        print(msg.format(path))
        sys.exit(1)

    s = Regex(r'SampleSheet_.+_R.+_L[0-9]+_*.*.csv')
    msg = 'ERROR: Samples Sheet file not named correctly.'
    msg += ' REQUIRED FORMAT: SampleSheet_{SEQUENCER}_R{RUN}_L{LANES}.csv'    
    Schema(s, error=msg).validate(path)

def write_errors(errors):
    """Writing out all errors
    """
    for section in errors.keys():
        if len(errors[section]) > 0:
            print('#-- Errors in section: {} --#'.format(section))
            x = [str(e) for  e in errors[section]]
            print('\n---\n'.join(x))
        else:
            continue
    
def main(args):
    """Main interface
    """
    validate_path(args.samples_sheet)
    sheet = read_samples_sheet(args.samples_sheet)
    errors = validate_samples_sheet(sheet, args.seq_tech)
    # errors?
    if any([len(x) > 0 for x in errors.values()]):
        write_errors(errors)
        sys.exit(1)
    else:
        print('Success! Sample sheet validation passed!')
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
