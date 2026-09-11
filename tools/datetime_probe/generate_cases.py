#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
from pathlib import Path
import re,subprocess,argparse
parser=argparse.ArgumentParser(description="Extract source DateTime literal and executable probe fixtures")
parser.add_argument('--source-root',type=Path,required=True)
parser.add_argument('--probe',type=Path,required=True)
parser.add_argument('--work-dir',type=Path,required=True)
parser.add_argument('--output-dir',type=Path,required=True)
args=parser.parse_args()
args.work_dir.mkdir(parents=True,exist_ok=True)
args.output_dir.mkdir(parents=True,exist_ok=True)
source=args.source_root/'src/tests/class_tests/openms/source/DateTime_test.cpp' 
lines=source.read_text().splitlines()
primary=[]
for i,line in enumerate(lines):
 m=re.search(r'(\w+)\.set\("([^"\\]*)"\)',line)
 if not m:continue
 var,value=m.groups()
 if 'TEST_EXCEPTION' in line:primary.append((i+1,value,'error',''));continue
 following='\n'.join(lines[i+1:i+4])
 expected=re.search(r'TEST_EQUAL\('+var+r'\.get\(\),\s*"([^"]*)"\)',following)
 if expected:primary.append((i+1,value,'ok',expected[1]))
(args.output_dir/'datetime_class_literals.tsv').write_text('source_line\tinput\tstatus\texpected_get\n'+''.join('\t'.join(map(str,row))+'\n' for row in primary))
cases=[]
def add(op,text='',fmt='',numbers='',seed=''):
 cases.append([f'case_{len(cases):04}',seed.encode().hex(),op,text.encode().hex(),fmt.encode().hex(),str(numbers)])
for _,value,_,_ in primary:add('set',value)
for value in ['', '2000-01-02','2000-01-02xZ','2000-01-02xZtail','2000-01-02 Z','2000-01-02Z ',
 '2000-01-02T03:04:05tail','2000-01-02T03:04:05-99:99','2000-01-02T03:04:05+bad',
 '+2000-01-02T03:04:05','2000-+1-02T03:04:05','2000-01-02T03:04:05junk.dot',
 '2000-01-02 03:04:05.1','2000-01-02T03:04:05/hi','2000-01-02\x00Z','2000-01-02Z\x00tail',
 '2000-01-02T03:04:05\x00tail','2000-01-02T03:04:05\x00.dot','2000-01-02T03:04:05\x00+',
 '2000-01-02T03:04:05\x00/', '2000-01-0203:04:05','2000-01-02 03:04:05 ignored',
 '0001-01-01T00:00:00','0000-01-01T00:00:00','2147483647-12-31T23:59:59',
 '1900-02-29T00:00:00','2000-02-29T00:00:00','2100-02-29T00:00:00',
 'Wed Dec 14 11:59:58 2006','Any Dec 14 2006','XXXDec142006','Any Dec 14 11:bad',
 'Any XXX 14 11:00:00 2006','Any December 14 2006','Any dec 14 2006','Any Dec 14 11:00',
 ' 2000-01-02T 03:\t04:\v05', '2000-01-02\v03:04:05']:
 add('set',value,seed='2011-08-05T15:32:07.468')
for digits in ['0','1','01','12','123','1234','12345','2147483647','-1',' -1','+1',' 1','-0']:
 for suffix in ['', '+02:00','Z']:
  add('set','2011-08-05T15:32:07.'+digits+suffix)
formats=['yyyy-MM-ddThh:mm:ss','yyyy-MM-ddThh:mm:ss.zzz','yyyy-MM-dd hh:mm:ss','yyyy-MM-dd+hh:mm','yyyy-MM-ddThh:mm:ssZ','yyyy-MM-dd','hh:mm:ss','unknown']
for fmt in formats:
 for value in ['2020-03-09T08:07:06','2020-03-09T08:07:06Z','2020-03-09T08:07:06.4',
 '2020-03-09T08:07:06.12345','2020-03-09T08:07:06. -1','2020-03-09T08:07:06.+1',
 '2020-03-09 08:07:06','2020-03-09+08:07','2020-03-09','08:07:06','24:00:00','2020-02-30']:
  add('from',value,fmt)
for op,inputs in [('date',['12/14/2006','14.12.2006','2006-12-14','2006-12-14ignored','1/2/2000-tail','2000-02-30','2000-1-1\x00-tail','-1.1.2000']),('time',['11:59:58','11:59:58tail','11:59:58\x00bad','+11: +59:58','24:00:00','23:60:00','23:59:60'])]:
 for value in inputs:
  add(op,value);add(op,value,seed='2011-08-05T15:32:07.468')
for seed in ['', 'time:11:59:58','0001-01-01T00:00:00','0100-03-01T00:00:00','0400-03-01T00:00:00',
 '1900-03-01T00:00:00','2000-03-01T00:00:00','2100-03-01T00:00:00',
 '2020-12-31T23:59:59.468','2147483647-01-01T00:00:00']:
 for seconds in [-2147483648,-86401,-1,0,1,86401,2147483647]:
  if seed.startswith('2147483647') and seconds>86401:continue
  add('add',numbers=seconds,seed=seed)
for op,numbers in [('parts','5,4,666,3,2,1'),('parts','12,14,2006,11,59,58'),('parts','2,30,2000,0,0,0'),('parts','1,1,2147483647,23,59,59'),('dateparts','12,14,2006'),('dateparts','2,30,2000'),('timeparts','11,59,58'),('timeparts','24,0,0')]:
 add(op,numbers=numbers,seed='2011-08-05T15:32:07.468')
add('none');add('clear',seed='2011-08-05T15:32:07.468')
(args.work_dir/'cases.tsv').write_text(''.join('\t'.join(c)+'\n' for c in cases))
with (args.work_dir/'cases.tsv').open('rb') as data:
 p=subprocess.run([str(args.probe)],stdin=data,capture_output=True,check=True)
(args.work_dir/'stderr.log').write_bytes(p.stderr)
assert not p.stderr
header='id\tseed_hex\toperation\tinput_hex\tformat_hex\tnumbers\tstatus\tvalid\tcomponents\tiso_hex\tmilliseconds_hex\tspace_hex\tplus_hex\tz_hex\tdate_hex\ttime_hex\tget_hex\tget_date_hex\tget_time_hex\n'
(args.output_dir/'datetime_cpp_probe.tsv').write_bytes(header.encode()+p.stdout)
print(f'{len(primary)} published class-test literal rows; {len(cases)} executed C++ oracle cases')
