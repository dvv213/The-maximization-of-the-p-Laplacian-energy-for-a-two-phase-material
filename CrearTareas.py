"""
Created on Wed Jan 20 16:42:42 2021

@author: donato
"""

import os 
import glob
import shutil
from MaximizationOfPLaplacianEnergy import *
print('Creando directorio soluciones')
save_dir='Tareas/'

if os.path.isdir(save_dir[:-1])==True:
    shutil.rmtree(save_dir[:-1])
os.mkdir(save_dir[:-1])
#shutil.copy2('ParallelVersion.py',save_dir[:-1] )   
#shutil.copy2('ParallelVersion2.py',save_dir[:-1] )

shutil.copy2('MaximizationOfPLaplacianEnergy.py',save_dir[:-1])
save_dir_sols='Soluciones/'
if os.path.isdir(save_dir_sols[:-1])==False:
    os.mkdir(save_dir_sols[:-1])
    ###Vaciamos la carpeta raiz


files=(glob.glob(save_dir_sols+'*'))
for file in files:
    os.remove(file)

c=1
k=1
def save_mesh(N):
    mesh=create_mesh(N)
    file_mesh=File('MallaN'+str(N)+'.xml')
    file_mesh<<mesh
def create_task(i,p,N,metodo):
    aux='import os\nprint(os.getcwd())'+\
'\nfrom MaximizationOfPLaplacianEnergy import *'+\
        '\nf = Expression("1",degree=0)'+\
    '\np='+str(p)+'\nN='+str(N)+'\nk='+str(k)+'\nc='+str(c)+\
'\nalpha=0.001*min(1-1/p,1/2)'+\
    '\nmetodo='+"'"+metodo+"'"+\
       "\nmesh=Mesh('MallaN'+str(N)+'.xml')"+\
       "\nprint('Solving problem for p',p,'N=',N)"+\
        '\nV,Theta= create_spaces(mesh)'+\
         '\nprint(Function(V).vector()[:].shape)'+\
        '\nsolution=solve_problem(Theta,V,p,c,k,f,metodo,max_it=2000,eps=10**-20,alpha=alpha,beta=0.5,max_it_init=100,tol=10**-7)'+\
        "\nimport dill\ndill.dump(solution,open('Soluciones/solution'+metodo.replace(' ','_')+'P'+str(p).replace('.','_')+'N'+str(N)+'.pkl','wb'))"
    with open(save_dir+'Tarea_'+str(i)+'.py','w') as file:
        file.write(aux)
#ps=[2]
#Ns=[40]
ps=[1.2,2,100]
Ns=[i for i in range(20,60,5)]
#ps=[1.1]
#Ns=[25]
metodos=['alternating optimization','feasible directions']
#ips=[10]
#Ns=[10]
for N in Ns:
    print('Creando Malla N:',N)
    save_mesh(N)
i=1
for p in ps:
    for N in Ns:
        for metodo in metodos:
            print('Creando Tarea p=',p,'N=',N,'metodo:',metodo)
            create_task(i,p,N,metodo)
            i+=1
NumeroTareas=len(ps)*len(Ns)*len(metodos)
print('Numero de Tareas:',NumeroTareas)
file=open('RefEjecutarTareas.sh'  ,'r' )
contenido=file.read().replace('NumeroTareas',str(NumeroTareas))
file.close()
file=open('EjecutarTareas.sh','w')
file.write(contenido)
file.close()
