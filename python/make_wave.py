
# +PATCH,//WAVE/PYTHON
# +DECK,wave_make,T=PYTHON.

import os
import sys
import platform
import glob

def Quit(*args, delay=0):
  #reakpoint()
  nargs =  len(args)

  text = ''
  for i in range(nargs):
    text += str(args[i]) + " "
  #endif

  if delay > 0:

    if len(text):
      print("\n",text, "\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    else:
      print("\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    #endif len(text):

    if platform.system() == 'Windows':
      os.system("sleep " + str(delay) + " && taskkill /F /PID " + str(os.getpid()) + " &")
    else:
      os.system("sleep " + str(delay) + " && kill " + str(os.getpid()) + " &")
    #endif platform.system() == 'Windows'

  elif delay < 0:
    return
  else:
    print("\n",text)
    if platform.system() == 'Windows':
      os.system("taskkill /F /PID " + str(os.getpid()))
    else:
      os.system("kill " + str(os.getpid()))
    #endif platform.system() == 'Windows'

#enddef Quit(text = '', delay=0)

def touch(filename):
  global OS, Iverbose, Idry
  path = os.getcwd()
  #reakpoint()
  if Iverbose: print('Touching ' + filename)
  if OS == 'Windows':
    scom = 'type ' + filename + ' > ' + path + '\\Kopie'
    if Idry == 0: os.system(scom)
    scom = 'type ' + path + '\\Kopie > ' + filename
    if Idry == 0: os.system(scom)
    scom = 'del ' + path + '\\Kopie'
    if Idry == 0: os.system(scom)
  else:
    scom = 'touch ' + filename
    if Idry == 0: os.system(scom)
  #endif
#enddef touch(file)

global Iverbose,Idry,Idebug,WI, Sepp, OS

args=sys.argv; nargs = len(args)

if platform.system() == 'Windows':
  Sepp = '\\'
  OS = 'Windows'
else:
  Sepp = '/' 
  OS = 'Linux'
#endif

WI = os.getcwd() + Sepp

Iverbose = 0
Idebug = 0
Idry = 0

if nargs > 1:
  try:
    Iverbose = int(args[1])
  except:
    n = '\n'
    print(n)
    print("Usage: python3 " + args[0] + " [verbose level]",n)
    print("To force total recompilation delete ",n,"bin" + Sepp + "wave.exe",n)
    Quit()
  #end try
#endif

if nargs > 2: Idebug = int(args[2])

global Wave_tree,Scomp_all,Scomp_omp,Scomp,Texe,Tlib,Scomp_nowarn

Scomp = "gfortran -std=legacy -c -O2 -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_nowarn = "gfortran -w -std=legacy -c -O2 -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_all = "gfortran -std=legacy -c -O2 -cpp -fcheck=all -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_omp = "gfortran -std=legacy -c -O2 -cpp -finit-local-zero -fcheck=all -fopenmp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -ffixed-line-length-none -funroll-loops "

try:
  import config_fortran as cf
  
  iconf = 0

  scomp = Scomp
  try:
    Scomp = cf.Ucomp
    iconf += 1
  except:
    Scomp = scomp
  #endtry
  
  scomp = Scomp_nowarn
  try:
    Scomp_nowarn = cf.Ucomp_nowarn
    iconf += 1
  except:
    Scomp_nowarn = scomp
  #endtry
  
  scomp = Scomp_all
  try:
    Scomp_all = cf.Ucomp_all
    iconf += 1
  except:
    Scomp_all = scomp
  #endtry
  
  scomp = Scomp_omp
  try:
    Scomp_omp = cf.Ucomp_omp
    iconf += 1
  except:
    Scomp_omp = scomp
  #endtry
  
  if iVerbose*iconf:
    print("\n FORTRAN command and options read taken from " + Sepp + "python" + Sepp + "config_fortran.py")
  #endif
except:
  pass
#endtry

def print_wave_tree():
  global Wave_tree

  ftree = open('wave.tre','w')
  itop = 0

  for topd in Wave_tree:

    itop+=1
    ftree.write('\n' + '---- Directory ' + str(itop) + ' ' + topd[0] + ' ' + str(topd[1]) + '\n\n')

    modfor = topd[2]
    imodfor = 0
    for mf in modfor:
      imodfor+=1
      slin = str(imodfor) + ' mod' + Sepp  + mf[0] + ' ' + str(mf[1])
      #print(slin)
      ftree.write(slin + '\n')
    #endfor

    ftree.write('\n')
    modmod = topd[3]
    imodmod = 0
    for mm in modmod:
      imodmod+=1
      slin = str(imodmod) + ' ' + mm[0] + ' ' + str(mm[1])
      #print(slin)
      ftree.write(slin + '\n')
    #endfor

    ftree.write('\n')
    cmn = topd[4]
    icmn = 0
    for cm in cmn:
      icmn+=1
      slin = str(icmn) + ' ' + cm[0] + ' ' + str(cm[1])
      #print(slin)
      ftree.write(slin + '\n')
    #endfor

    ftree.write('\n')
    ff = topd[5]
    iff = 0
    for f in ff:
      iff+=1
      slin = str(iff) + ' ' + f[0] + ' ' + str(f[1])
      #print(slin)
      ftree.write(slin + '\n')
    #endfor

  #endfor

  ftree.close()
#    Wave_tree.append([topd,t,modfor,modmod,cmn,fort])
#enddef print_wave_tree()

def get_wave_tree():

  global WI,Wave_tree,Iverbose,Idry,Idebug,Texe,Tlib,Sepp

  try:
    Texe = os.stat(WI + Sepp + 'bin' + Sepp + 'wave.exe').st_mtime_ns
  except:
    Texe = 0
  #endtry

  top = glob.glob(WI+"*")

  Wave_tree = []
  #reakpoint()

  for topd in top:

    dd = topd.split(Sepp)[-1]

    if dd == 'cmz' or dd == 'doc' or dd == 'check_system' or dd == 'bin' \
    or dd == 'python' or dd == 'main' or dd == 'lib': continue

    t = os.stat(topd).st_mtime_ns

    modf = glob.glob(topd+Sepp+"mod"+Sepp+"*.f")

    modfor = []
    for ff in modf:
      f = ff.split(Sepp)[-1]
      tf = os.stat(ff).st_mtime_ns
      modfor.append([f,tf])
    #endfor

    modm = glob.glob(topd+Sepp+"*.mod")
    modmod = []
    for ff in modm:
      f = ff.split(Sepp)[-1]
      tf = os.stat(ff).st_mtime_ns
      modmod.append([f,tf])
    #endfor

    cm = glob.glob(topd+Sepp+"*.cmn")
    cmn = []
    for ff in cm:
      f = ff.split(Sepp)[-1]
      tf = os.stat(ff).st_mtime_ns
      cmn.append([f,tf])
    #endfor

    ff = glob.glob(topd+Sepp+"*.f")
    #reakpoint()
    fort = []
    for fff in ff:
      f = fff.split(Sepp)[-1]
      tf = os.stat(fff).st_mtime_ns
      fort.append([f,tf])
    #endfor

    Wave_tree.append([topd,t,modfor,modmod,cmn,fort])

  #endfor get_wave_tree

  #print_wave_tree()

#enddef get_wave_tree

def wave_update():

  global WI,Wave_tree,Texe,Scomp_all,Scomp_omp,Scomp,Iverbose,Idry,Idebug,Scomp_nowarn,Sepp

  kmain = 0

  #reakpoint()
  get_wave_tree()
  #reakpoint()
  for td in Wave_tree:

    dd = td[0]
    ds = dd + Sepp
    dsm = dd + Sepp + "mod" + Sepp
    t = td[1]
    modfor = td[2]
    cmn = td[4]
    #reakpoint()
    fort = td[5]

    scomp = Scomp

    lib = ''
    libm = ''
    ranl = 0
    ranlm = 0
    slibm = []
    slib = []

    ddd = dd.split(Sepp)[-1]

    #reakpoint()

    if ddd == 'mhbook':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libmhbook.a'
      libm = WI + 'lib' + Sepp + 'libmhbook_modules.a'
    elif ddd == 'mshcern':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libmshcern.a'
      libm = WI + 'lib' + Sepp + 'libmshcern_modules.a'
      scomp = Scomp_nowarn
    elif ddd == 'mshplt':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libmshplt.a'
      libm = WI + 'lib' + Sepp + 'libmshplt_modules.a'
    elif ddd == 'nomp':
      #reakpoint()
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libwave.a'
      libm = WI + 'lib' + Sepp + 'libwave_modules.a'
      scomp = Scomp_all
    elif ddd == 'omp':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libwave_omp.a'
      libm = WI + 'lib' + Sepp + 'libwave_omp_modules.a'
      scomp = Scomp_omp
    elif ddd == 'urad':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'liburad.a'
      libm = WI + 'lib' + Sepp + 'liburad_modules.a'
      scomp = Scomp_all  # uradcfft does boundary tricks
    elif ddd == 'user':
      if Iverbose >= 0: print("\nProcessing",dd)
      lib = WI + 'lib' + Sepp + 'libuser.a'
      libm = WI + 'lib' + Sepp + 'libuser_modules.a'
      scomp = Scomp_omp
    #endif
    #reakpoint()

    klib = 0
    Tlib = Texe + 1
    try:
      Tlib = os.stat(lib).st_mtime_ns
      if Tlib > Texe: klib = 1
    except:
      klib = 1
    #endtry

    klibm = 0
    Tlib = Texe + 1
    try:
      Tlib = os.stat(libm).st_mtime_ns
      if Tlib > Texe: klibm = 1
    except:
      klibm = 1
    #endtry

    scompmod = "cd " + dd + Sepp + "mod && " + scomp + '-J.. '
    scomp = "cd " + dd + " && " + scomp
    
    for f in modfor: # Compile modules

      ff = f[0]
      t = f[1]

      if t < Texe and klibm == 0: continue

      if Iverbose > 0: print(ff)

      fo = ff[:-1] + "o"
      fm = ff[:-1] + "mod"

      Flines = open(ds+"mod"+Sepp+ff,'r')

      while True:
        l = Flines.readline()
        if not l: break
        sl = l.split()
        if len(sl) == 0: continue
        key = sl[0].lower()
        if key== 'module':
          m = sl[1].lower()
          break
        #endif
      #end while
      Flines.close()

      if Iverbose > 0: print("\nModule:",m)
      scom = scompmod + "-o " + fo + " " + ff
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      
      #reakpoint()
#      if platform.system() == 'Windows':
#        scom = 'move ' + dsm + m + ".mod " + dd
#      else:
#        scom = 'mv ' + dsm + m + ".mod " + dd
#      #endif
      
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)

      slibm.append(dsm + fo)
      ranlm = 1

      # Search use of module in *.cmn
      for ft in cmn:

        f = ft[0]
        t = ft[1]

        if t < Texe and klibm == 0: continue

        Flines = open(ds+f,'r')
        if Idebug > 1: print("\n",ds+f)

        while True:
          l = Flines.readline()
          if Idebug > 1: print(l)
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          if len(sl) > 1:
            key = sl[0].lower()
            if key== 'use':
              if sl[1].lower() == m:
                # scom = 'touch ' + ds+f
                #if Iverbose > 0: print("\n",scom,"\n")
                #if Idry == 0: os.system(scom)
                scom = touch(ds+f)
                break
              #endif
            #endif
          #endif
        #end while
        Flines.close()

      #endfor

      # Search use of module in *.f
      for ft in fort:

        f = ft[0]
        t = ft[1]

        if Idebug > 1: print(f)
        #reakpoint()
        Flines = open(ds+f,'r')
        while True:
          l = Flines.readline()
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          if len(sl) > 1:
            key = sl[0].lower()
            if key== 'implicit':
              if sl[1].lower() == 'none': break
            elif key== 'use':
              if sl[1].lower() == m:
                #scom = 'touch ' + ds+f
                #if Iverbose > 0: print("\n",scom,"\n")
                #if Idry == 0: os.system(scom)
                scom = touch(ds+f)
                break
              #endif
            #endif
          #endif
        #end while
        Flines.close()

      #endfor

    #endfor modfor

    if ranlm:
      #reakpoint()
      for ob in slibm:
        scom = 'ar rc ' + libm + " " + ob
        if Idry == 0: os.system(scom)
      #endif
      
      scom = 'ranlib ' + libm
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      ranlm = 0
      slibm = []
      kmain = 1
    #endif

    # Check *.cmn

    for ft in cmn:

      f = ft[0]
      #print(f)
      #if f == 'genfun.cmn': debug()
      t = os.stat(ds+f).st_mtime_ns

      if t < Texe and klib == 0: continue

      fcmn = f.split(Sepp)[-1]

      for fft in fort:

        #if fft[0] == 'erzfun.f': debug('erzfun')

        Flines = open(ds+fft[0],'r')
        #Flines = open(ds+fft[0],'r',errors='ignore')
        #Flines = open(ds+fft[0],'r',encoding='latin1') #, errors='ignore')
#        nlin = 0
        while True:
          l = Flines.readline()
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          #if fft[0] == 'erzfun.f':
#            nlin += 1
#            print(nlin,l)
#            if nlin > 80: debug(nlin)
          if len(sl) < 2: continue
          #print(sl)
          if sl[0][0] == '*' or sl[0][0] == '!' or len(sl[0]) < 7: continue
          key = sl[0].lower()
          #if fft[0] == 'erzfun.f': print(sl)
          if key== 'include':
#            if fft[0] == 'erzfun.f': debug('include')
            if sl[1].lower() == "'" + fcmn + "'" or sl[1].lower() == '"' + fcmn + '"':
              #scom = 'touch ' + ds+fft[0]
              #if Iverbose > 0: print("\n",scom,"\n")
              #if Idry == 0: os.system(scom)
              scom = touch(ds+f)
              break
            #endif
          #endif
        #end while
        Flines.close()
#        if fft[0] == 'erzfun.f': Quit()
      #endfor fort

    #endfor cmn

    # Compile *.f if neccessary

    for ft in fort:

      f = ft[0]
#      print(f)
      t = os.stat(ds+f).st_mtime_ns

      if t < Texe and klib == 0: continue

      fo = f[:-1] + "o"

      scom = scomp + "-o " + fo + " " + f
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)

      slib.append(ds + fo)
      ranl = 1

    #endfor

    if ranl:
      for ob in slib:
        scom = 'ar rc ' + lib + " " + ob
        if Idry == 0: os.system(scom)
      #endif
#      scom = 'ar rc ' + lib + " " + slib
#      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      scom = 'ranlib ' + lib
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      kmain = 1
      ranl = 0
      slib = []
    #endif

  #endfor dir

  if kmain:
    #scom = WI + Sepp + "shell" + Sepp + "compile_wave_incl.sh"
    #if Iverbose > 0: print("\n",scom,"\n")
    #if Idry == 0: os.system(scom)
    #if Iverbose >=0: print("\n--- " + WI  + "bin" + Sepp + "wave.exe updated ---\n")
    if Idry == 0: wave_compile()
  else:
    if Iverbose >=0: print("\n--- No need to update " + WI  + "bin" + Sepp + "wave.exe ---\n")
  #endif


#enddef wave_update

def debug(key='debug'):
  print("debug:",key)
#endif

def wave_compile():
  
  global WI,Sepp,Iverbose,Idry,Move,Delete
  
  #Idry = 1
  
  if Iverbose:
    print("\nMaking .." + Sepp + "bin" + Sepp + "wave.exe\n")

  if platform.system() == 'Windows':
    Move = 'move '
    Delete = 'del '
  else:
    Move = 'mv '
    Delete = 'rm '
  #endif
  
  pathmain = 'main' + Sepp
  pathmod = pathmain + 'mod' + Sepp
  
  #scom = Delete + pathmod + '*.o'

#  if Idry: print("\n")
#  if Idry: print(scom,"\n")
#  else: os.system(scom)

  #scom = Delete + pathmod + '*.mod'
  
  #if Idry: print(scom,"\n")
  #else: os.system(scom)

  scom = 'cd ' + pathmod + ' && ' \
  'gfortran -c -O2' + \
  ' -fcheck=bounds -fbacktrace' + \
  ' -ffpe-summary=invalid,zero,overflow' + \
  ' -fdec -fd-lines-as-comments' + \
  ' -Wno-align-commons -fno-automatic -ffixed-line-length-none' + \
  ' -finit-local-zero -J..' + \
  ' -funroll-loops *.f'
  
  if Idry: print(scom,"\n")
  else: os.system(scom)

  #scom = move + pathmod + '*.mod ' + pathmain
  
  #if Idry: print(scom,"\n")
  #else: os.system(scom)

  scom = 'cd main && ' \
  'gfortran -O2 -cpp' + \
  ' -fopenmp' + \
  ' -fcheck=bounds' + \
  ' -fbacktrace' + \
  ' -ffpe-summary=invalid,zero,overflow' + \
  ' -fdec -fd-lines-as-comments' + \
  ' -Wno-align-commons' + \
  ' -ffixed-line-length-none' + \
  ' -finit-local-zero' + \
  ' -funroll-loops' + \
  ' -o ..' + Sepp + 'bin' + Sepp + 'wave.exe wave_main.f' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libuser.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libuser_modules.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libmhbook.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libmhbook_modules.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libmshplt.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libmshcern.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave_modules.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave_omp.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave_omp_modules.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave.a' + \
  ' ..' + Sepp + 'lib' + Sepp + 'libwave_omp.a'
  #endif 

  if Idry: print(scom,"\n")
  else: os.system(scom)

#enddef wave_compile()  

wave_update()
