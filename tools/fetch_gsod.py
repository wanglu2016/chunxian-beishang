import os, re, sys, shutil, subprocess, py7zr, json
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
OUT=os.path.join(BASE,'parsed'); os.makedirs(OUT,exist_ok=True); os.makedirs(os.path.join(BASE,'raw'),exist_ok=True); os.makedirs(os.path.join(BASE,'x'),exist_ok=True)
pat=re.compile(r'^\d{4}/(45|46|5\d)\d{4}-')
y0,y1=int(sys.argv[1]),int(sys.argv[2])
for y in range(y0,y1+1):
    outf=os.path.join(OUT,f'{y}.tsv')
    if os.path.exists(outf): continue
    arc=os.path.join(BASE,'raw',f'{y}.7z')
    if not os.path.exists(arc):
        r=subprocess.run(['curl','-s','-f','-m','600','-o',arc,f'https://raw.githubusercontent.com/CalebBell/gsod/master/{y}.7z'])
        if r.returncode!=0: print('FAIL download',y,flush=True); continue
    xd=os.path.join(BASE,'x',str(y)+'_tmp'); shutil.rmtree(xd,ignore_errors=True); os.makedirs(xd)
    with py7zr.SevenZipFile(arc) as z: names=[n for n in z.getnames() if pat.match(n)]
    with py7zr.SevenZipFile(arc) as z: z.extract(path=xd,targets=names)
    rows=[]
    for root,_,fs in os.walk(xd):
        for f in fs:
            if not f.endswith('.op'): continue
            for line in open(os.path.join(root,f),errors='replace'):
                t=line.split()
                if len(t)<18 or not t[0].isdigit(): continue
                stn,wban,d,temp,cnt=t[0],t[1],t[2],t[3],t[4]
                mx,mn=t[17].rstrip('*'),t[18].rstrip('*')
                rows.append(f'{stn}\t{wban}\t{d}\t{temp}\t{cnt}\t{mx}\t{mn}\n')
    with open(outf+'.part','w') as o: o.writelines(rows)
    os.rename(outf+'.part',outf)
    shutil.rmtree(xd); os.remove(arc)
    print('done',y,len(names),'files',len(rows),'rows',flush=True)
