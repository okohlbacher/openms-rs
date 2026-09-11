# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
"""Independent direct translation of SDK82 IntegerMassDecomposer table loops.
No Rust import/output and no C++ execution. For bounded small positive inputs only.
"""
from math import gcd
from itertools import combinations

def table(weights):
    first,second=weights[:2]; infinity=first*weights[-1]
    rows=[[infinity]*first for _ in weights]
    for row in rows: row[0]=0
    witnesses=[(0,0)]*first
    increment=second%first; residue=increment; mass=second; count=0
    while residue:
        rows[1][residue]=mass; mass+=second; count+=1
        witnesses[residue]=(1,count)
        residue=(residue+increment)%first
    for i,current_mass in enumerate(weights[2:],2):
        if current_mass>=rows[i-1][current_mass%first]:
            rows[i]=rows[i-1].copy(); continue
        d=gcd(first,current_mass); previous=rows[i-1]; current=rows[i]
        if d==1:
            value=residue=count=0
            for _ in range(first):
                value+=current_mass; residue=(residue+current_mass)%first; count+=1
                if value>previous[residue]: value=previous[residue];count=0
                else: witnesses[residue]=(i,count)
                current[residue]=value
        else:
            cur=current_mass%first; prev=0; increment=cur-d; counters=[0]*first
            current[1:d]=previous[1:d]
            for _ in range(first//d-1):
                for _ in range(d):
                    counters[cur]+=1
                    if current[prev]+current_mass>previous[cur]:current[cur]=previous[cur];counters[cur]=0
                    else:current[cur]=current[prev]+current_mass;witnesses[cur]=(i,counters[cur])
                    prev+=1;cur+=1
                prev=cur-d;cur=(cur+increment)%first
            cont=True; loops=0
            while cont:
                loops+=1;assert loops<1000
                cont=False;prev+=1;cur+=1;counters[cur]+=1
                for _ in range(1,d):
                    if current[prev]+current_mass<current[cur]:
                        current[cur]=current[prev]+current_mass;cont=True
                        witnesses[cur]=(i,counters[cur])
                    else:counters[cur]=0
                    prev+=1;cur+=1
                prev=cur-d;cur=(cur+increment)%first
    return infinity,rows,witnesses

if __name__=='__main__':
    import json
    records=[]
    for weights,mass in [([10,6,15],33),([10,16,25],73)]:
        infinity,rows,witnesses=table(weights)
        r=mass%weights[0]
        brute=[]
        def visit(i,left,result):
            if i==len(weights)-1:
                if left%weights[i]==0:brute.append(result+[left//weights[i]])
                return
            for n in range(left//weights[i]+1):visit(i+1,left-n*weights[i],result+[n])
        visit(0,mass,[])
        records.append({'weights':weights,'mass':mass,'infinity':infinity,'rows':rows,
          'queried_witness':witnesses[r],'exists':rows[-1][r]!=infinity and mass>=rows[-1][r],
          'mathematical_compositions':brute})
    assert records[1]['rows'][-1] == [0,41,32,73,64,25,16,57,48,89]
    assert records[1]['queried_witness'] == (2,0)
    assert records[1]['exists'] and records[1]['mathematical_compositions'] == [[0,3,1]]
    # C++ getDecomposition: m=73, r=3, i=2, j=0; subtracting j*25
    # leaves m and r unchanged, so while(m!=0) cannot terminate.
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[2]/'tests/data/ims_witness_source_oracle.json')
    args=parser.parse_args()
    text=json.dumps(records,indent=2)+'\n'
    if args.check:
        if args.output.read_text()!=text:raise SystemExit('reference differs from source recurrence')
        print('Verified two source-recurrence rows and independent exact compositions.')
    else:args.output.write_text(text)
