"""Independent source-expression oracle; not a Rust or full SDK execution."""
from pathlib import Path
import itertools,json,math,re,struct,sys
if len(sys.argv) != 2:
 raise SystemExit('usage: generate_metabo_peptide_reference.py PINNED_SDK')
root=Path(sys.argv[1])
p=root/'src/openms/source/CHEMISTRY/ElementDB.cpp';source=p.read_text()
f32=lambda x:struct.unpack('f',struct.pack('f',x))[0]
symbols={'C':'carbon','H':'hydrogen','N':'nitrogen','O':'oxygen','S':'sulfur','P':'phosphorus'}
maps={}
for sym,name in symbols.items():
 maps[sym]={}
 for key in ['abundance','mass']:
  text=re.search(r'\b'+name+'_'+key+r' = (.*);',source)[1]
  maps[sym][key]={int(i):float(x) for i,x in re.findall(r'\{(\d+)u,\s*([^}]+)\}',text)}
avg={s:sum(m['mass'][i]*v for i,v in m['abundance'].items()) for s,m in maps.items()}
comp={'C':4.9384,'H':7.7583,'N':1.3577,'O':1.4773,'S':.0417,'P':0.}
factor=1000/sum(comp[s]*avg[s] for s in 'CHNOSP')
counts={s:math.floor(comp[s]*factor+.5) for s in 'CNOSP'}
counts['H']=math.floor((1000-sum(counts[s]*avg[s] for s in 'CNOSP'))/avg['H']+.5)

def convolve(a,b,limit,cast):
 n=min(len(a)+len(b)-1,limit)
 out=[0.]*n
 for i in range(len(a)-1,-1,-1):
  for j in range(min(n-i,len(b))-1,-1,-1):
   out[i+j]=cast(out[i+j]+cast(a[i]*b[j]))
 return out

def power(a,n,limit,cast):
 if n==1:return a[:]
 out=a[:] if n&1 else [1.]
 square=convolve(a,a,limit+1,cast)
 i=1
 while True:
  if n&(1<<i):out=convolve(out,square,limit,cast)
  if i>=(n-1).bit_length():break
  square=convolve(square,square,limit+1,cast)
  i+=1
 return out

def envelope(order,limit,cast):
 out=[1.]
 for s in order:
  ab=maps[s]['abundance'];a=[cast(ab.get(i,0)) for i in range(min(ab),max(ab)+1)]
  out=convolve(out,power(a,counts[s],limit,cast),limit,cast)
 total=0.
 for v in reversed(out):total+=v
 return [cast(v/total) for v in out]

def score(a,h):
 x=[v/max(a) for v in a];y=[v/max(h) for v in h]
 xy=xx=yy=0.
 for u,v in zip(x,y):xy+=u*v;xx+=u*u;yy+=v*v
 return xy/(math.sqrt(xx)*math.sqrt(yy))
rows=[]
for h in ([10.,10.],[10.,5.],[10.,10.,10.],[10.,5.,2.]):
 actual=[envelope(order,len(h),f32) for order in itertools.permutations('HCNO')]
 direct=envelope('HCNO',len(h),f32)
 native=envelope('CHNO',len(h),float)
 rows.append({'hypothesis':h,'source_f32_HCNO_distribution':direct,'source_f32_HCNO_score':score(direct,h),'source_f32_all24orders_distribution_min':[min(r[i] for r in actual) for i in range(len(h))],'source_f32_all24orders_distribution_max':[max(r[i] for r in actual) for i in range(len(h))],'source_f32_all24orders_score_interval':[min(score(r,h) for r in actual),max(score(r,h) for r in actual)],'independent_f64_distribution':native,'independent_f64_score':score(native,h)})
result={'origin':'Independent Python stdlib translation of source expressions, not source class literal or C++/Rust output. All24 active element permutations account for source pointer-key formula traversal. Zero S/P and charge0 identity convolutions do not change values.','weight':1000.,'formula_counts':counts,'factor':factor,'rows':rows}
print(json.dumps(result,indent=2))
