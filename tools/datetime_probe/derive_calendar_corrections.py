#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Independent month stepping using Python calendar.monthrange, no Rust code/output.

Avoids the native March-era day-number algorithm. Only the explicitly identified
macOS timegm early-input-year cases are projected. All raw C++ rows stay intact.
"""
import calendar
import argparse
from pathlib import Path
parser=argparse.ArgumentParser(description="Derive platform-independent DateTime arithmetic references by calendar month stepping")
parser.add_argument('--data-dir',type=Path,required=True)
args=parser.parse_args()
rows=[]
def move_month(y,m,amount):
 m+=amount
 if m==0:return y-1,12
 if m==13:return y+1,1
 return y,m
for line in (args.data_dir/'datetime_cpp_probe.tsv').read_text().splitlines()[1:]:
 f=line.split('\t');seed=bytes.fromhex(f[1]).decode()
 if f[2]!='add' or seed not in ['', 'time:11:59:58','0001-01-01T00:00:00','0100-03-01T00:00:00','0400-03-01T00:00:00']:continue
 if seed=='':y,m,d,h,minute,s=0,0,0,0,0,0
 elif seed.startswith('time:'):y,m,d,h,minute,s=0,0,0,11,59,58
 else:
  date,time=seed.split('T');y,m,d=map(int,date.split('-'));h,minute,s=map(int,time.split(':'))
 if m==0:y-=1;m=12
 if d==0:y,m=move_month(y,m,-1);d=calendar.monthrange(y,m)[1]
 days,remainder=divmod(h*3600+minute*60+s+int(f[5]),86400)
 d+=days
 while d<1:
  y,m=move_month(y,m,-1);d+=calendar.monthrange(y,m)[1]
 while d>calendar.monthrange(y,m)[1]:
  d-=calendar.monthrange(y,m)[1];y,m=move_month(y,m,1)
 h,remainder=divmod(remainder,3600);minute,s=divmod(remainder,60)
 date=f'{y:04}-{m:02}-{d:02}';time=f'{h:02}:{minute:02}:{s:02}'
 valid=f[7]=='1'
 formats=[date+'T'+time,date+'T'+time+'.000',date+' '+time,date+f'+{h:02}:{minute:02}',date+'T'+time+'Z',date,time]
 get=[date+' '+time,date,time]
 if not valid:formats=['']*7;get=['0000-00-00 00:00:00','0000-00-00','00:00:00']
 rows.append([f[0],','.join(map(str,[m,d,y,h,minute,s]))]+[v.encode().hex() for v in formats+get])
assert len(rows)==35
(args.data_dir/'datetime_calendar_corrections.tsv').write_text('id\tcomponents\tiso_hex\tmilliseconds_hex\tspace_hex\tplus_hex\tz_hex\tdate_hex\ttime_hex\tget_hex\tget_date_hex\tget_time_hex\n'+''.join('\t'.join(r)+'\n' for r in rows))
print('35 independent Python calendar month-step correction rows written')
