import hashlib,json,pathlib,sqlite3,zlib
root=pathlib.Path('/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R')
base=root/'.reference/openms4-core-bc9cc12/src/tests/class_tests/openms/data'
p=base/'SqliteMassFile_1.sqMass'
c=sqlite3.connect(p.resolve().as_uri()+'?mode=ro',uri=True)
out={'evidence':'Python sqlite3 read-only fixture inspection and independently executed SQL, not executed C++','sqlite_version':sqlite3.sqlite_version,'fixture':{'path':str(p.relative_to(root)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size},'tables':{},'indexes':list(c.execute("select name,tbl_name,sql from sqlite_master where type='index' and sql is not null"))}
for name in ['RUN','SPECTRUM','CHROMATOGRAM','PRECURSOR','PRODUCT']:out['tables'][name]=list(c.execute('select * from '+name))
out['arrays']=list(c.execute('select SPECTRUM_ID,CHROMATOGRAM_ID,COMPRESSION,DATA_TYPE,length(DATA) from DATA'))
out['run_extra']=[{'run_id':r[0],'compressed_bytes':len(r[1]),'xml_bytes':len(zlib.decompress(r[1]))} for r in c.execute('select RUN_ID,DATA from RUN_EXTRA')]
c.close()
c=sqlite3.connect(':memory:')
c.executescript('CREATE TABLE PRECURSOR(SPECTRUM_ID INT, CHROMATOGRAM_ID INT, ISOLATION_TARGET REAL,ISOLATION_LOWER REAL,ISOLATION_UPPER REAL); CREATE TABLE SPECTRUM(ID INT PRIMARY KEY, NATIVE_ID TEXT,MSLEVEL INT); CREATE TABLE DATA(SPECTRUM_ID INT); INSERT INTO SPECTRUM VALUES (0,\'s0\',2),(1,\'s1\',2); INSERT INTO PRECURSOR VALUES(NULL,0,412.5,12.5,12.5),(0,NULL,412.5,12.5,12.5),(1,NULL,412.5,10,10); INSERT INTO DATA VALUES(1),(0);')
out['swath_null_row']=list(c.execute('SELECT SPECTRUM_ID FROM PRECURSOR WHERE ISOLATION_TARGET BETWEEN 412.49 AND 412.51'))
out['swath_distinct_tuples']=list(c.execute('SELECT DISTINCT(ISOLATION_TARGET),ISOLATION_TARGET-ISOLATION_LOWER,ISOLATION_TARGET+ISOLATION_UPPER FROM PRECURSOR INNER JOIN SPECTRUM ON SPECTRUM_ID=SPECTRUM.ID WHERE MSLEVEL==2'))
out['metadata_order']=list(c.execute('SELECT SPECTRUM.ID,SPECTRUM.NATIVE_ID FROM SPECTRUM'))
out['blob_join_order']=list(c.execute('SELECT SPECTRUM.ID,SPECTRUM.NATIVE_ID FROM SPECTRUM INNER JOIN DATA ON SPECTRUM.ID=DATA.SPECTRUM_ID'))
print(json.dumps(out,indent=2))
