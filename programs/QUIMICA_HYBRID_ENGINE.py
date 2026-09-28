# Motor MicroPython de QUIMICA_HYBRID 2.8.0
# Pegar este contenido en la pestaña Python si fuera necesario.
import hpprime as hp

# Valores didacticos tomados del HPP de referencia del usuario.
H=6.63e-34
C=3.0e8
EV=1.6e-19
ERG=1.0e7
A0=0.529
V0=2.18e8
RY=13.6
SYMBOLS=('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
         'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
         'Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr '
         'Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()

def gv(name):
    return hp.eval(name)

def load_constants():
    global H,C,EV,ERG,A0,V0,RY
    pairs=(('CFG_H',6.63e-34),('CFG_C',3.0e8),('CFG_EV',1.6e-19),
           ('CFG_ERG',1.0e7),('CFG_A0',0.529),('CFG_V0',2.18e8),('CFG_RY',13.6))
    vals=[]
    for name,default in pairs:
        try:
            v=float(gv(name)); vals.append(v if v>0 else default)
        except: vals.append(default)
    H,C,EV,ERG,A0,V0,RY=vals

load_constants()

def sv(name,value):
    # El parser HOME de algunos firmwares falla con exponentes e+NN.
    num=repr(float(value)).replace('e+','e')
    hp.eval(name+':='+num)

def st(text):
    # Construye los saltos con CHAR(10), compatible con PPL.
    parts=str(text).split('\n')
    quoted=[]
    for part in parts:
        part=part.replace('\\','\\\\').replace('"','\\"')
        quoted.append('"'+part+'"')
    hp.eval('QTEXT:='+'+CHAR(10)+'.join(quoted))

def ok(a=0,b=0,c=0,text=''):
    hp.eval('QERR:=0')
    sv('QOUT1',a); sv('QOUT2',b); sv('QOUT3',c)
    st(text)

def fail(text):
    hp.eval('QERR:=1')
    st(text)

def sci(x):
    if x==0: return '0'
    s='{:.6g}'.format(x)
    if 'e' in s:
        a,b=s.split('e'); n=int(b); digs='⁰¹²³⁴⁵⁶⁷⁸⁹'
        exp=('⁻' if n<0 else '')+''.join(digs[int(d)] for d in str(abs(n)))
        return a+'×10'+exp
    return s

def planck(op,x):
    if x<=0: raise ValueError('El dato debe ser mayor que cero.')
    if op==101: lam=x*1e-9
    elif op==102: lam=C/x
    else: lam=H*C/x
    f=C/lam; e=H*f
    return f,lam,e

def bohr(op,n1,n2):
    if n1<1 or int(n1)!=n1: raise ValueError('n debe ser entero positivo.')
    n1=int(n1)
    if op==201:
        return A0*n1*n1,V0/n1,-RY/(n1*n1)
    if n2<1 or int(n2)!=n2 or n1==n2:
        raise ValueError('Use dos niveles enteros positivos y distintos.')
    n2=int(n2)
    de=RY*(1/(n1*n1)-1/(n2*n2))
    ej=abs(de)*EV
    return de,ej,H*C/ej

def einstein(op,x):
    if x<0: raise ValueError('El dato no puede ser negativo.')
    if op==301:
        kg=x; e=kg*C*C
    elif op==302:
        kg=x/1000.0; e=kg*C*C
    else:
        e=x; kg=e/(C*C)
    return kg,e,e/EV

def atom(z,a,charge):
    if z<1 or a<z or int(z)!=z or int(a)!=a or int(charge)!=charge:
        raise ValueError('Compruebe Z, A y la carga entera.')
    p=int(z); n=int(a-z); e=int(z-charge)
    if e<0: raise ValueError('La carga produce electrones negativos.')
    return p,n,e

def electron_config(z):
    if z<1 or z>118 or int(z)!=z:
        raise ValueError('Z debe ser entero entre 1 y 118.')
    order=((1,'s',2),(2,'s',2),(2,'p',6),(3,'s',2),(3,'p',6),(4,'s',2),
           (3,'d',10),(4,'p',6),(5,'s',2),(4,'d',10),(5,'p',6),(6,'s',2),
           (4,'f',14),(5,'d',10),(6,'p',6),(7,'s',2),(5,'f',14),(6,'d',10),(7,'p',6))
    left=int(z); filled=[]
    for n,l,cap in order:
        if left<=0: break
        e=min(left,cap); filled.append([n,l,e]); left-=e
    # Excepciones de estabilidad mas usadas en cursos introductorios.
    shifts={24:(3,'d',4,'s'),29:(3,'d',4,'s'),41:(4,'d',5,'s'),42:(4,'d',5,'s'),
            44:(4,'d',5,'s'),45:(4,'d',5,'s'),46:(4,'d',5,'s'),47:(4,'d',5,'s'),
            78:(5,'d',6,'s'),79:(5,'d',6,'s')}
    if z in shifts:
        dn,dl,sn,sl=shifts[z]
        di=-1; si=-1
        for indice in range(len(filled)):
            x=filled[indice]
            if x[0]==dn and x[1]==dl: di=indice
            if x[0]==sn and x[1]==sl: si=indice
        if di>=0 and si>=0 and filled[si][2]>0:
            filled[di][2]+=1
            filled[si][2]-=1
    active=[]
    for x in filled:
        if x[2]>0: active.append(x)
    partes=[]
    supers='⁰¹²³⁴⁵⁶⁷⁸⁹'
    for x in active:
        exponente=''
        for digito in str(x[2]):
            exponente+=supers[int(digito)]
        partes.append(str(x[0])+x[1]+exponente)
    cfg=' '.join(partes)
    period=0
    for x in filled:
        if x[2]>0 and x[0]>period: period=x[0]
    last_n,last_l,last_e=active[-1]
    if last_l=='s': group=last_e
    elif last_l=='p': group=12+last_e
    elif last_l=='d': group=2+last_e
    else: group=3
    if z==2: group=18
    valence=0
    for x in active:
        if x[0]==period: valence+=x[2]
    if 57<=z<=71: family='Lantanido'; block='f'
    elif 89<=z<=103: family='Actinido'; block='f'
    elif group==1: family='Alcalino' if z!=1 else 'No metal'
    elif group==2: family='Alcalinoterreo'
    elif group==17: family='Halogeno'
    elif group==18: family='Gas noble'
    elif 3<=group<=12: family='Metal de transicion'
    else: family='Familia del grupo '+str(group)
    if not (57<=z<=71 or 89<=z<=103): block=last_l
    return SYMBOLS[int(z)-1],cfg,period,group,valence,family,block

def gcd(a,b):
    a=abs(int(a)); b=abs(int(b))
    while b: a,b=b,a%b
    return a or 1

def roman(n):
    vals=((10,'X'),(9,'IX'),(5,'V'),(4,'IV'),(1,'I'))
    out=''
    for v,s in vals:
        while n>=v: out+=s; n-=v
    return out

def prefix(n):
    p=('','mono','di','tri','tetra','penta','hexa','hepta','octa','nona','deca')
    return p[n] if n<len(p) else str(n)

# Estados positivos usados por la nomenclatura tradicional escolar. Los
# elementos sin serie tradicional se nombran de forma funcional, sin inventar
# terminaciones. Esto permite aceptar los 118 elementos con una salida honesta.
OX_STATES={
 'H':(1,), 'Li':(1,), 'Be':(2,), 'B':(3,), 'C':(2,4), 'N':(1,2,3,4,5),
 'Na':(1,), 'Mg':(2,), 'Al':(3,), 'Si':(2,4), 'P':(1,3,5), 'S':(2,4,6),
 'Cl':(1,3,5,7), 'K':(1,), 'Ca':(2,), 'Sc':(3,), 'Ti':(2,3,4),
 'V':(2,3,4,5), 'Cr':(2,3,6), 'Mn':(2,3,4,6,7), 'Fe':(2,3),
 'Co':(2,3), 'Ni':(2,3), 'Cu':(1,2), 'Zn':(2,), 'Ga':(3,), 'Ge':(2,4),
 'As':(3,5), 'Se':(2,4,6), 'Br':(1,3,5,7), 'Rb':(1,), 'Sr':(2,),
 'Y':(3,), 'Zr':(4,), 'Nb':(3,5), 'Mo':(2,3,4,5,6), 'Tc':(4,7),
 'Ru':(2,3,4,6,8), 'Rh':(1,3), 'Pd':(2,4), 'Ag':(1,), 'Cd':(2,),
 'In':(1,3), 'Sn':(2,4), 'Sb':(3,5), 'Te':(2,4,6), 'I':(1,3,5,7),
 'Cs':(1,), 'Ba':(2,), 'La':(3,), 'Ce':(3,4), 'Pr':(3,4), 'Nd':(3,),
 'Sm':(2,3), 'Eu':(2,3), 'Gd':(3,), 'Tb':(3,4), 'Dy':(3,), 'Ho':(3,),
 'Er':(3,), 'Tm':(2,3), 'Yb':(2,3), 'Lu':(3,), 'Hf':(4,), 'Ta':(5,),
 'W':(2,3,4,5,6), 'Re':(2,3,4,5,6,7), 'Os':(2,3,4,6,8),
 'Ir':(1,3,4), 'Pt':(2,4), 'Au':(1,3), 'Hg':(1,2), 'Tl':(1,3),
 'Pb':(2,4), 'Bi':(3,5), 'Po':(2,4,6), 'Fr':(1,), 'Ra':(2,),
 'Ac':(3,), 'Th':(4,), 'Pa':(4,5), 'U':(3,4,5,6), 'Np':(3,4,5,6,7),
 'Pu':(3,4,5,6), 'Am':(2,3,4,5,6), 'Cm':(3,4), 'Bk':(3,4),
 'Cf':(2,3,4), 'Es':(2,3), 'Fm':(2,3), 'Md':(2,3), 'No':(2,3), 'Lr':(3,),
 'Rf':(4,), 'Db':(5,), 'Sg':(6,), 'Bh':(7,), 'Hs':(8,), 'Mt':(1,3),
 'Ds':(2,4), 'Rg':(1,3), 'Cn':(2,), 'Nh':(1,3), 'Fl':(2,4),
 'Mc':(1,3), 'Lv':(2,4), 'Ts':(1,3,5,7), 'Og':(2,)}

TR_ROOT={
 'H':'hidr','Li':'lit','Be':'beril','B':'bor','C':'carbon','N':'nitr',
 'Na':'sod','Mg':'magnes','Al':'alumin','Si':'silic','P':'fosfor','S':'sulfur',
 'Cl':'clor','K':'potas','Ca':'calc','Sc':'escand','Ti':'titan','V':'vanad',
 'Cr':'crom','Mn':'mangan','Fe':'ferr','Co':'cobalt','Ni':'niquel','Cu':'cupr',
 'Zn':'cinc','Ga':'gal','Ge':'german','As':'arsen','Se':'selen','Br':'brom',
 'Rb':'rubid','Sr':'estronc','Y':'itr','Zr':'circon','Nb':'niob','Mo':'molibd',
 'Tc':'tecnec','Ru':'ruten','Rh':'rod','Pd':'palad','Ag':'argent','Cd':'cadm',
 'In':'ind','Sn':'estan','Sb':'antimon','Te':'telur','I':'yod','Cs':'ces',
 'Ba':'bar','La':'lantan','Ce':'cer','Pr':'praseodim','Nd':'neodim','Sm':'samar',
 'Eu':'europ','Gd':'gadolin','Tb':'terb','Dy':'dispros','Ho':'holm','Er':'erb',
 'Tm':'tul','Yb':'iterb','Lu':'lutec','Hf':'hafn','Ta':'tantal','W':'wolfram',
 'Re':'ren','Os':'osm','Ir':'irid','Pt':'platin','Au':'aur','Hg':'mercur',
 'Tl':'tal','Pb':'plumb','Bi':'bismut','Po':'polon','Fr':'franc','Ra':'rad',
 'Ac':'actin','Th':'tor','Pa':'protactin','U':'uran','Np':'neptun','Pu':'pluton',
 'Am':'americ','Cm':'cur','Bk':'berkel','Cf':'californ','Es':'einsten',
 'Fm':'ferm','Md':'mendelev','No':'nobel','Lr':'laurenc'}

NONMETALS=set(('B C N Si P S Cl As Se Br Te I').split())
SPECIAL_HYDRIDES={'B':'borano','C':'metano','N':'amoniaco','Si':'silano',
                  'P':'fosfina','As':'arsina','Sb':'estibina'}

def trad_adj(sym,val):
    states=OX_STATES.get(sym,())
    root=TR_ROOT.get(sym,'')
    if not states or val not in states or not root: return ''
    if len(states)==1: return root+'ico'
    i=states.index(val); n=len(states)
    if n==2: return root+('oso' if i==0 else 'ico')
    if n==3: return ('hipo'+root+'oso' if i==0 else root+'oso' if i==1 else root+'ico')
    # En series de cuatro o mas, la tradicional escolar solo fija con claridad
    # los dos extremos y los estados intermedios bajo/alto.
    if i==0: return 'hipo'+root+'oso'
    if i==1: return root+'oso'
    if i==n-2: return root+'ico'
    if i==n-1: return 'per'+root+'ico'
    return ''

def traditional(sym,val,kind,name='',anion=''):
    if kind=='hidracido': return anion
    if kind=='hidruro' and sym in SPECIAL_HYDRIDES: return SPECIAL_HYDRIDES[sym]
    if kind=='oxido' and sym=='H': return 'agua'
    adj=trad_adj(sym,val)
    if not adj:
        return 'Sin nombre tradicional estandar; use Stock: '+name+' ('+roman(val)+')'
    if kind=='oxido': return ('anhidrido ' if sym in NONMETALS else 'oxido ')+adj
    if kind=='peroxido': return 'peroxido '+adj
    if kind=='hidruro': return 'hidruro '+adj
    if kind=='hidroxido': return 'hidroxido '+adj
    if kind=='binario': return anion+' '+adj
    if kind=='oxosal': return anion+' '+adj
    return adj

RADICALS={
  1:('OH',-1,'hidroxido',{'O':1,'H':1}), 2:('NO3',-1,'nitrato',{'N':1,'O':3}),
  3:('NO2',-1,'nitrito',{'N':1,'O':2}), 4:('SO4',-2,'sulfato',{'S':1,'O':4}),
  5:('SO3',-2,'sulfito',{'S':1,'O':3}), 6:('CO3',-2,'carbonato',{'C':1,'O':3}),
  7:('PO4',-3,'fosfato',{'P':1,'O':4}), 8:('ClO',-1,'hipoclorito',{'Cl':1,'O':1}),
  9:('ClO2',-1,'clorito',{'Cl':1,'O':2}), 10:('ClO3',-1,'clorato',{'Cl':1,'O':3}),
 11:('ClO4',-1,'perclorato',{'Cl':1,'O':4}), 12:('MnO4',-1,'permanganato',{'Mn':1,'O':4}),
 13:('CrO4',-2,'cromato',{'Cr':1,'O':4}), 14:('Cr2O7',-2,'dicromato',{'Cr':2,'O':7}),
 15:('CN',-1,'cianuro',{'C':1,'N':1}), 16:('NH4',1,'amonio',{'N':1,'H':4}),
 17:('BrO',-1,'hipobromito',{'Br':1,'O':1}), 18:('BrO2',-1,'bromito',{'Br':1,'O':2}),
 19:('BrO3',-1,'bromato',{'Br':1,'O':3}), 20:('BrO4',-1,'perbromato',{'Br':1,'O':4}),
 21:('IO',-1,'hipoyodito',{'I':1,'O':1}), 22:('IO2',-1,'yodito',{'I':1,'O':2}),
 23:('IO3',-1,'yodato',{'I':1,'O':3}), 24:('IO4',-1,'peryodato',{'I':1,'O':4}),
 25:('BO3',-3,'borato',{'B':1,'O':3}), 26:('SiO3',-2,'silicato',{'Si':1,'O':3}),
 27:('AsO3',-3,'arsenito',{'As':1,'O':3}), 28:('AsO4',-3,'arsenato',{'As':1,'O':4}),
 29:('SeO3',-2,'selenito',{'Se':1,'O':3}), 30:('SeO4',-2,'selenato',{'Se':1,'O':4}),
 31:('MnO4',-2,'manganato',{'Mn':1,'O':4}), 32:('S2O3',-2,'tiosulfato',{'S':2,'O':3}),
 33:('C2O4',-2,'oxalato',{'C':2,'O':4}), 34:('CH3COO',-1,'acetato',{'C':2,'H':3,'O':2}),
 35:('HCO3',-1,'hidrogenocarbonato (bicarbonato)',{'H':1,'C':1,'O':3}),
 36:('HSO4',-1,'hidrogenosulfato (bisulfato)',{'H':1,'S':1,'O':4}),
 37:('PO3',-3,'fosfito',{'P':1,'O':3}), 38:('HPO4',-2,'hidrogenofosfato',{'H':1,'P':1,'O':4}),
 39:('H2PO4',-1,'dihidrogenofosfato',{'H':2,'P':1,'O':4}),
 40:('HSO3',-1,'hidrogenosulfito (bisulfito)',{'H':1,'S':1,'O':3}),
 41:('OCN',-1,'cianato',{'O':1,'C':1,'N':1}), 42:('SCN',-1,'tiocianato',{'S':1,'C':1,'N':1}),
 43:('O2',-2,'peroxido',{'O':2}), 44:('O2',-1,'superoxido',{'O':2}),
 45:('S2O8',-2,'peroxodisulfato (persulfato)',{'S':2,'O':8}),
 46:('TeO3',-2,'telurito',{'Te':1,'O':3}), 47:('TeO4',-2,'telurato',{'Te':1,'O':4}),
 48:('MoO4',-2,'molibdato',{'Mo':1,'O':4}), 49:('WO4',-2,'wolframato',{'W':1,'O':4}),
 50:('VO3',-1,'metavanadato',{'V':1,'O':3}), 51:('VO4',-3,'ortovanadato',{'V':1,'O':4}),
 52:('AlO2',-1,'aluminato',{'Al':1,'O':2}), 53:('H3O',1,'hidronio',{'H':3,'O':1}),
 54:('HCOO',-1,'formiato',{'H':1,'C':1,'O':2}),
 55:('C6H5O7',-3,'citrato',{'C':6,'H':5,'O':7})}

BINARY_ANIONS={'F':'fluoruro','Cl':'cloruro','Br':'bromuro','I':'yoduro','At':'astaturo',
 'S':'sulfuro','Se':'seleniuro','Te':'telururo','N':'nitruro','P':'fosfuro',
 'As':'arseniuro','Sb':'antimoniuro','C':'carburo','Si':'siliciuro','B':'boruro',
 'H':'hidruro','O':'oxido'}

HYDRACIDS={'F':'acido fluorhidrico','Cl':'acido clorhidrico','Br':'acido bromhidrico',
 'I':'acido yodhidrico','At':'acido astatohidrico','S':'acido sulfhidrico',
 'Se':'acido selenhidrico','Te':'acido telurhidrico','CN':'acido cianhidrico'}

OXOACIDS={'NO3':'acido nitrico','NO2':'acido nitroso','SO4':'acido sulfurico',
 'SO3':'acido sulfuroso','CO3':'acido carbonico','PO4':'acido fosforico',
 'ClO':'acido hipocloroso','ClO2':'acido cloroso','ClO3':'acido clorico',
 'ClO4':'acido perclorico','BrO':'acido hipobromoso','BrO2':'acido bromoso',
 'BrO3':'acido bromico','BrO4':'acido perbromico','IO':'acido hipoyodoso',
 'IO2':'acido yodoso','IO3':'acido yodico','IO4':'acido peryodico',
 'BO3':'acido borico','SiO3':'acido metasilicico','AsO3':'acido arsenioso',
 'AsO4':'acido arsenico','SeO3':'acido selenioso','SeO4':'acido selenico',
 'MnO4':'acido manganico o permanganico segun la carga','CrO4':'acido cromico',
 'Cr2O7':'acido dicromico','S2O3':'acido tiosulfurico','C2O4':'acido oxalico'}
OXOACIDS.update({'PO3':'acido fosforoso','TeO3':'acido teluroso','TeO4':'acido telurico',
 'MoO4':'acido molibdico','WO4':'acido tungstico','VO3':'acido metavanadico',
 'VO4':'acido ortovanadico','OCN':'acido cianico','SCN':'acido tiocianico',
 'S2O8':'acido peroxodisulfurico'})

def counts_text(parts):
    order=[]; total={}
    for comp,mul in parts:
        for s,n in comp.items():
            if s not in total: order.append(s); total[s]=0
            total[s]+=n*mul
    return ', '.join(s+'='+str(total[s]) for s in order)

def add_count(dst,src,mul=1):
    for k,v in src.items(): dst[k]=dst.get(k,0)+v*mul

def parse_formula(raw):
    raw=str(raw); charge=0
    if raw.endswith('}') and '{' in raw:
        k=raw.rfind('{'); charge=int(raw[k+1:-1]); raw=raw[:k]
    s=raw; stack=[{}]; i=0
    while i<len(s):
        ch=s[i]
        if ch in '([': stack.append({}); i+=1
        elif ch in ')]':
            if len(stack)<2: raise ValueError('Parentesis incorrectos en '+raw)
            group=stack.pop(); i+=1; j=i
            while j<len(s) and s[j].isdigit(): j+=1
            mul=int(s[i:j] or '1'); add_count(stack[-1],group,mul); i=j
        elif ch.isupper():
            j=i+1
            while j<len(s) and s[j].islower(): j+=1
            el=s[i:j]
            if el not in SYMBOLS: raise ValueError('Simbolo desconocido: '+el)
            k=j
            while k<len(s) and s[k].isdigit(): k+=1
            stack[-1][el]=stack[-1].get(el,0)+int(s[j:k] or '1'); i=k
        else: raise ValueError('Caracter no valido en '+raw)
    if len(stack)!=1: raise ValueError('Falta cerrar parentesis en '+raw)
    return raw,stack[0],charge

def fgcd(a,b): return gcd(a,b)
def frac(n,d=1):
    if d==0: raise ValueError('Division por cero.')
    if d<0: n=-n; d=-d
    g=fgcd(n,d); return (int(n)//g,int(d)//g)
def fadd(a,b): return frac(a[0]*b[1]+b[0]*a[1],a[1]*b[1])
def fneg(a): return (-a[0],a[1])
def fmul(a,b): return frac(a[0]*b[0],a[1]*b[1])
def fdiv(a,b): return frac(a[0]*b[1],a[1]*b[0])
def lcm(a,b): return abs(a*b)//fgcd(a,b)

def null_vector(mat):
    if not mat or not mat[0]: raise ValueError('Ecuacion vacia.')
    a=[[(int(x),1) for x in row] for row in mat]
    rows=len(a); cols=len(a[0]); piv=[]; r=0
    for c in range(cols):
        p=next((i for i in range(r,rows) if a[i][c][0]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]; lead=a[r][c]
        a[r]=[fdiv(x,lead) for x in a[r]]
        for i in range(rows):
            if i!=r and a[i][c][0]:
                q=a[i][c]; a[i]=[fadd(a[i][j],fneg(fmul(q,a[r][j]))) for j in range(cols)]
        piv.append(c); r+=1
        if r==rows: break
    free=[j for j in range(cols) if j not in piv]
    if not free: raise ValueError('Solo existe la solucion trivial; revise la ecuacion.')
    x=[(0,1) for _ in range(cols)]
    for j in free: x[j]=(1,1)
    for i in range(len(piv)-1,-1,-1):
        c=piv[i]; total=(0,1)
        for j in free: total=fadd(total,fmul(a[i][j],x[j]))
        x[c]=fneg(total)
    den=1
    for q in x: den=lcm(den,q[1])
    out=[q[0]*(den//q[1]) for q in x]
    if all(v<0 for v in out): out=[-v for v in out]
    if any(v<=0 for v in out):
        raise ValueError('La ecuacion es ambigua o requiere agregar especies del medio.')
    g=out[0]
    for v in out[1:]: g=fgcd(g,v)
    return [v//g for v in out]

SUBS='₀₁₂₃₄₅₆₇₈₉'
SUPS='⁰¹²³⁴⁵⁶⁷⁸⁹'
def sub_number(n): return ''.join(SUBS[int(x)] for x in str(n))
def sup_number(n):
    sign='⁺' if n>0 else '⁻'
    return sign+''.join(SUPS[int(x)] for x in str(abs(n)))
def oxnum(n): return ('+'+str(n)) if n>0 else ('−'+str(abs(n)) if n<0 else '0')
def pretty_formula(f):
    out=''; i=0
    while i<len(f):
        if f[i].isdigit():
            j=i
            while j<len(f) and f[j].isdigit(): j+=1
            out+=sub_number(f[i:j]); i=j
        else: out+=f[i]; i+=1
    return out
def show_species(item):
    f,c,q=item; pf=pretty_formula(f)
    return pf if q==0 else '['+pf+']'+sup_number(q)

KNOWN_NAMES={'KNO2':'nitrito de potasio','K2O':'oxido de potasio / oxido potasico',
 'O2':'oxigeno molecular / dioxigeno','NO':'oxido de nitrogeno (II) / oxido nitrico',
 'Zn':'zinc','NO3':'nitrato','OH':'hidroxido','H2O':'agua','Zn(OH)4':'tetrahidroxozincato(II)',
 'NH3':'amoniaco','FeSO4':'sulfato de hierro(II) / sulfato ferroso',
 'H2SO4':'acido sulfurico','HNO3':'acido nitrico','Fe2(SO4)3':'sulfato de hierro(III) / sulfato ferrico',
 'HCl':'cloruro de hidrogeno / acido clorhidrico','Cl2':'cloro molecular'}

def oxidation_states(item):
    f,comp,q=item
    if len(comp)==1: return {next(iter(comp)):0}
    st={}; fixed={'H':1,'O':-2,'F':-1,'Li':1,'Na':1,'K':1,'Rb':1,'Cs':1,
      'Be':2,'Mg':2,'Ca':2,'Sr':2,'Ba':2,'Al':3}
    for e in comp:
        if e in fixed: st[e]=fixed[e]
    patterns=(('SO4','S',6),('SO3','S',4),('NO3','N',5),('NO2','N',3),
      ('CO3','C',4),('PO4','P',5),('PO3','P',3),('ClO4','Cl',7),
      ('ClO3','Cl',5),('ClO2','Cl',3),('ClO','Cl',1))
    for pat,e,v in patterns:
        if pat in f and e in comp: st[e]=v
    unknown=[e for e in comp if e not in st]
    if len(unknown)==1:
        e=unknown[0]; known=sum(comp[k]*v for k,v in st.items())
        num=q-known
        if num%comp[e]==0: st[e]=num//comp[e]
    return st

def balance_equation(encoded):
    if '>' not in encoded: raise ValueError('Falta la flecha de reaccion.')
    ls,rs=encoded.split('>',1)
    left=[parse_formula(x) for x in ls.split('|') if x]
    right=[parse_formula(x) for x in rs.split('|') if x]
    if not left or not right: raise ValueError('Complete reactivos y productos.')
    elems=[]
    for it in left+right:
        for e in it[1]:
            if e not in elems: elems.append(e)
    mat=[]
    for e in elems:
        mat.append([x[1].get(e,0) for x in left]+[-x[1].get(e,0) for x in right])
    mat.append([x[2] for x in left]+[-x[2] for x in right])
    coef=null_vector(mat); nl=len(left)
    def side(items,cs):
        return ' + '.join(('' if c==1 else str(c))+show_species(x) for x,c in zip(items,cs))
    original=side(left,[1]*nl)+' → '+side(right,[1]*len(right))
    final=side(left,coef[:nl])+' → '+side(right,coef[nl:])
    # Cambios de oxidacion que el motor puede deducir sin ambiguedad.
    lvals={}; rvals={}
    for x in left:
        for e,v in oxidation_states(x).items(): lvals.setdefault(e,[]).append(v)
    for x in right:
        for e,v in oxidation_states(x).items(): rvals.setdefault(e,[]).append(v)
    changes=[]
    for e in elems:
        av=sorted(set(lvals.get(e,[]))); bv=sorted(set(rvals.get(e,[])))
        pairs=[]
        for va in av:
            for vb in bv:
                if va!=vb and (va not in bv or vb not in av): pairs.append((va,vb))
        for va,vb in pairs:
            typ='OXIDACIÓN' if vb>va else 'REDUCCIÓN'
            line=e+': '+oxnum(va)+' → '+oxnum(vb)+'  '+typ+' ('+str(abs(vb-va))+' e⁻ por átomo)'
            if line not in changes: changes.append(line)
    change_text='\n'.join(changes) if changes else 'No se detectó un cambio redox único; se aplicó balance algebraico.'
    checks=[]
    for e in elems:
        a=sum(c*x[1].get(e,0) for x,c in zip(left,coef[:nl])); b=sum(c*x[1].get(e,0) for x,c in zip(right,coef[nl:]))
        checks.append(e+': '+str(a)+' = '+str(b))
    ql=sum(c*x[2] for x,c in zip(left,coef[:nl])); qr=sum(c*x[2] for x,c in zip(right,coef[nl:]))
    nomen=[]
    for x in left+right: nomen.append(show_species(x)+': '+KNOWN_NAMES.get(x[0],'nombre disponible en el modulo de nomenclatura'))
    return ('ECUACIÓN INGRESADA\n'+original+'\n\n1. ASIGNACIÓN AUTOMÁTICA\nEl motor aplica H=+1, O=−2, elementos libres=0 y la suma igual a la carga.\n\n2. OXIDACIÓN Y REDUCCIÓN\n'+change_text+'\n\n3. BALANCE DE MATERIA Y CARGA\nSe construyó y resolvió el sistema de conservación para '+str(len(elems))+' elementos.\n\n4. ECUACIÓN BALANCEADA\n'+final+'\n\n5. COMPROBACIÓN\n'+'\n'.join(checks)+'\nCarga: '+str(ql)+' = '+str(qr)+'\n\n6. NOMENCLATURA\n'+'\n'.join(nomen))

def group_compound(z,val,rid,mode):
    z=int(z); val=int(val); rid=int(rid); mode=int(mode)
    if rid not in RADICALS: raise ValueError('Radical no reconocido.')
    rf,rv,rname,rcomp=RADICALS[rid]
    names=hp.eval('TP_NOMBRES()')
    if mode==2:
        if rv>=0 or rf=='OH' or rid in (34,35,36,38,39,40,43,44,52,54,55):
            raise ValueError('Ese ion no se usa como oxoacido directo en este constructor.')
        s1='H'; name1='hidrogeno'; val=1; z=1
    else:
        if z<1 or z>118: raise ValueError('Seleccione un elemento valido.')
        s1=SYMBOLS[z-1]; name1=str(names[z-1])
    if val==0 or val*rv>=0: raise ValueError('Las cargas deben tener signos opuestos.')
    g=gcd(val,rv); n1=abs(rv)//g; ng=abs(val)//g
    grouped=len(rcomp)>1 and ng>1
    reverse=rv>0
    if reverse:
        formula=('('+rf+')' if grouped else rf)+('' if ng==1 else str(ng))+s1+('' if n1==1 else str(n1))
    else:
        formula=s1+('' if n1==1 else str(n1))+('('+rf+')' if grouped else rf)+('' if ng==1 else str(ng))
    q1=n1*val; q2=ng*rv
    if mode==2:
        kind='OXOACIDO'
        stock=('acido permanganico' if rid==12 else 'acido manganico' if rid==31 else OXOACIDS.get(rf,'acido de '+rname))
        systematic='hidrogeno + '+rname
        reason='Comienza con H y contiene un oxoanion.'
        tradname=stock
        if rf=='CN':
            kind='HIDRACIDO'; stock='cianuro de hidrogeno'
            systematic='cianuro de hidrogeno'; tradname='acido cianhidrico'
            reason='En disolucion acuosa se nombra con la terminacion -hidrico.'
    else:
        kind='HIDROXIDO' if rf=='OH' else ('SAL DE AMONIO' if rf=='NH4' else 'OXOSAL')
        if reverse:
            anname=BINARY_ANIONS.get(s1,name1.lower()+'uro')
            stock=anname+' de '+rname; systematic=stock; tradname=stock
            reason='Combina el cation '+rname+' con el anion '+anname+'.'
        else:
            stock=rname+' de '+name1+(' ('+roman(abs(val))+')' if s1!='H' else '')
            systematic=(prefix(ng)+rname+' de '+prefix(n1)+name1.lower()).replace('monohidroxido','hidroxido')
            reason=('Contiene el grupo OH unido a un cation.' if rf=='OH' else 'Combina un cation con un oxoanion.')
            tradname=traditional(s1,abs(val),'hidroxido' if rf=='OH' else 'oxosal',name1,rname)
    st_var('QSYM1',s1); st_var('QSYM2',rf); st_var('QGROUP',rf); st_var('QFORM',formula)
    sv('QSUB1',n1); sv('QSUB2',ng); sv('QV1',val); sv('QV2',rv)
    sv('QREV',1 if reverse else 0)
    atoms=counts_text([({s1:1},n1),(rcomp,ng)])
    return ('1. INTERPRETACIÓN\nFórmula: '+pretty_formula(formula)+'\nFunción: '+kind+'\n'+reason+'\n\n2. IONES O GRUPOS\n'+s1+': '+str(val)+'; '+pretty_formula(rf)+': '+str(rv)+'\n\n3. CRUCE Y SIMPLIFICACIÓN\nMCD = '+str(g)+'\nSubíndices: '+str(n1)+' y '+str(ng)+('\nEl paréntesis conserva unido al radical.' if grouped else '')+'\n\n4. COMPROBACIÓN\n'+str(n1)+'('+str(val)+') + '+str(ng)+'('+str(rv)+') = '+str(q1+q2)+'\nCarga total = 0\nÁtomos: '+atoms+'\n\n5. NOMENCLATURA\nStock: '+stock+'\nSistemática: '+systematic+'\nTradicional: '+tradname)

def formulate(z1,z2,v1,v2):
    z1=int(z1); z2=int(z2); v1=int(v1); v2=int(v2)
    if z1<1 or z1>118 or z2<1 or z2>118 or z1==z2:
        raise ValueError('Seleccione dos elementos diferentes.')
    if v1==0 or v2==0 or v1*v2>0:
        raise ValueError('Las valencias deben tener signos opuestos.')
    g=gcd(v1,v2); n1=abs(v2)//g; n2=abs(v1)//g
    s1=SYMBOLS[z1-1]; s2=SYMBOLS[z2-1]
    names=hp.eval('TP_NOMBRES()')
    name1=str(names[z1-1]); name2=str(names[z2-1])
    formula=s1+('' if n1==1 else str(n1))+s2+('' if n2==1 else str(n2))
    anion=''
    if s2=='O' and v2==-1:
        kind='peroxido'
        stock='peroxido de '+name1
        systematic=(prefix(n2)+'oxido de '+prefix(n1)+name1.lower()).replace('monooxido','monoxido')
    elif s2=='O' and v2==-2:
        kind='oxido'
        stock='oxido de '+name1+' ('+roman(abs(v1))+')'
        systematic=(prefix(n2)+'oxido de '+prefix(n1)+name1.lower()).replace('monooxido','monoxido')
    elif s2=='H' and v2==-1:
        kind='hidruro'
        stock='hidruro de '+name1+' ('+roman(abs(v1))+')'
        systematic=prefix(n2)+'hidruro de '+prefix(n1)+name1.lower()
    elif s1=='H' and v1==1:
        kind='hidracido'
        anion_name=BINARY_ANIONS.get(s2,name2.lower()+'uro')
        stock=anion_name+' de hidrogeno'
        systematic=prefix(n2)+anion_name+' de '+prefix(n1)+'hidrogeno'
        anion=HYDRACIDS.get(s2,'acido '+name2.lower()+'hidrico')
    else:
        kind='binario'
        anion=BINARY_ANIONS.get(s2,name2.lower()+'uro')
        stock=anion+' de '+name1+' ('+roman(abs(v1))+')'
        systematic=prefix(n2)+anion+' de '+prefix(n1)+name1.lower()
    trad=traditional(s1,abs(v1),kind,name1,anion)
    sv('QSUB1',n1); sv('QSUB2',n2); sv('QV1',v1); sv('QV2',v2)
    st_var('QSYM1',s1); st_var('QSYM2',s2)
    check=str(n1)+'('+str(v1)+') + '+str(n2)+'('+str(v2)+') = 0'
    text=('1. ELEMENTOS\n'+s1+' = '+name1+', oxidacion '+str(v1)+'\n'+
          s2+' = '+name2+', oxidacion '+str(v2)+'\n\n2. CRUCE DE VALENCIAS\n'
          +str(abs(v1))+' y '+str(abs(v2))+' se simplifican por '+str(g)+'\n'
          +'Subíndices: '+str(n1)+' y '+str(n2)+'\n\n3. FÓRMULA\n'+pretty_formula(formula)+
          '\n\n4. NOMENCLATURA\nStock: '+stock+'\nSistematica: '+systematic+
          '\nTradicional: '+trad+'\n\n5. COMPROBACION\n'+check)
    return text

def st_var(name,text):
    text=str(text).replace('\\','\\\\').replace('"','\\"')
    hp.eval(name+':="'+text+'"')

try:
    op=int(gv('QOP')); a=float(gv('QA')); b=float(gv('QB'))
    c=float(gv('QC')); d=float(gv('QD'))
    if 101<=op<=103:
        f,l,e=planck(op,a)
        if op==101:
            data='λ = '+sci(a)+' nm\nλ = '+sci(l)+' m'
        elif op==102:
            data='f = '+sci(a)+' Hz'
        else:
            data='E = '+sci(a)+' J'
        text=('1. DATOS\n'+data+'\n\n2. FORMULAS\n'
              'c = λ·ν\nE = h·ν\n\n3. SUSTITUCIÓN\n'
              'ν = c/λ\nν = '+sci(C)+' / '+sci(l)+'\n'
              'ν = '+sci(f)+' Hz\nE = ('+sci(H)+')('+sci(f)+')\n'
              'E = '+sci(e)+' J\n\n4. RESULTADO\n'
              'λ = '+sci(l)+' m\nν = '+sci(f)+' Hz\n'
              'E = '+sci(e)+' J\nE = '+sci(e/EV)+' eV')
        ok(f,l,e,text)
    elif 201<=op<=202:
        x,y,z=bohr(op,a,b)
        if op==201:
            text=('1. DATOS\nn = '+str(int(a))+'\n\n2. FORMULAS\n'
                  'rₙ = 0.529·n²\nvₙ = 2.18×10⁸/n\nEₙ = −13.6/n²\n\n'
                  '3. SUSTITUCIÓN\nrₙ = 0.529('+str(int(a))+')²\n'
                  'vₙ = 2.18×10⁸/'+str(int(a))+'\n'
                  'Eₙ = −13.6/'+str(int(a))+'²\n\n4. RESULTADO\n'
                  'Radio = '+sci(x)+' Å\nVelocidad = '+sci(y)+' cm·s⁻¹\nEnergía = '+sci(z)+' eV')
            ok(x,y,z,text)
        else:
            kind='EMISION' if x<0 else 'ABSORCION'
            text=('1. DATOS\nni = '+str(int(a))+'\nnf = '+str(int(b))+'\n\n'
                  '2. FÓRMULA\nΔE = 13.6(1/nᵢ² − 1/n_f²)\n'
                  'E(fotón) = |ΔE|\nλ = h·c/E\n\n3. SUSTITUCIÓN\n'
                  'ΔE = 13.6(1/'+str(int(a))+'² − 1/'+str(int(b))+'²)\n'
                  'ΔE = '+sci(x)+' eV\nE(fotón) = '+sci(y)+' J\n\n'
                  '4. RESULTADO\nTipo = '+kind+'\nλ = '+sci(z)+' m\n'
                  'λ = '+sci(z*1e9)+' nm\nν = '+sci(C/z)+' Hz')
            ok(x,y,z,text)
    elif 301<=op<=303:
        kg,e,ev=einstein(op,a)
        known=('m = '+sci(a)+' kg' if op==301 else
               'm = '+sci(a)+' g = '+sci(kg)+' kg' if op==302 else
               'E = '+sci(a)+' J')
        formula='E = m·c²' if op<303 else 'm = E/c²'
        subst=('E = ('+sci(kg)+')('+sci(C)+')²' if op<303 else
               'm = '+sci(e)+' / ('+sci(C)+')²')
        text=('1. DATOS\n'+known+'\nc = '+sci(C)+' m·s⁻¹\n\n2. FÓRMULA\n'+formula+
              '\n\n3. SUSTITUCIÓN\n'+subst+'\n\n4. RESULTADO\n'
              'm = '+sci(kg)+' kg\nm = '+sci(kg*1000)+' g\n'
              'E = '+sci(e)+' J\nE = '+sci(e*1e7)+' erg')
        ok(kg,e,ev,text)
    elif op==401:
        p,n,e=atom(a,b,c)
        species='NEUTRO' if c==0 else ('CATION' if c>0 else 'ANION')
        text=('1. DATOS\nZ = '+str(int(a))+'\nA = '+str(int(b))+'\nCarga = '+str(int(c))+
              '\n\n2. RELACIONES\np⁺ = Z\nn⁰ = A − Z\ne⁻ = Z − carga\n\n'
              '3. SUSTITUCIÓN\nn = '+str(int(b))+' − '+str(int(a))+'\n'
              'e⁻ = '+str(int(a))+' − ('+str(int(c))+')\n\n4. RESULTADO\n'
              'Protones = '+str(p)+'\nNeutrones = '+str(n)+'\nElectrones = '+str(e)+'\nEspecie = '+species)
        ok(p,n,e,text)
    elif op==402:
        sym,cfg,period,group,valence,family,block=electron_config(a)
        text=('1. ELEMENTO\n'+sym+'   Z = '+str(int(a))+'\n\n2. REGLA DE LLENADO\n'
              'Aufbau + Pauli + Hund\n\n3. CONFIGURACION\n'+cfg+'\n\n'
              '4. ANALISIS\nGrupo = '+str(group)+'\nPeriodo = '+str(period)+
              '\nElectrones de valencia = '+str(valence)+'\nFamilia = '+family+'\nBloque = '+block)
        ok(a,period,group,text)
    elif op==1001:
        text=formulate(gv('QZ1'),gv('QZ2'),gv('QV1'),gv('QV2'))
        ok(gv('QSUB1'),gv('QSUB2'),0,text)
    elif op==1002:
        text=group_compound(gv('QZ1'),gv('QV1'),gv('QRID'),gv('QMODE'))
        ok(gv('QSUB1'),gv('QSUB2'),0,text)
    elif op==1003:
        rid=int(gv('QRID'))
        if rid not in RADICALS: raise ValueError('Radical no reconocido.')
        st_var('QFORM',RADICALS[rid][0]); ok(RADICALS[rid][1],0,0,RADICALS[rid][2])
    elif op==2001:
        ok(0,0,0,balance_equation(gv('QRXN')))
    else: fail('Operacion del motor no reconocida.')
except Exception as ex:
    fail('Error: '+str(ex))
