
"""
Created on Wed Jan 20 17:19:37 2021

@author: donato
"""
from matplotlib import pyplot as plt
import os
import glob 
import dill
import numpy as np
from MaximizationOfPLaplacianEnergy import *
save_dir='Plots/'
if os.path.isdir(save_dir[:-1])==False:
    os.mkdir(save_dir[:-1])
    ###Vaciamos la carpeta raiz
    
files=(glob.glob(save_dir+'*'))
for file in files:
    os.remove(file)
    
save_dir_sols='Soluciones/'
files=(glob.glob(save_dir_sols+'sol*'))
solutions={}
ps=[]
Ns=[]
metodos=[]
for file in files:
    file_name=str(file)
    metodo=file_name[file_name.find('solution')+len('solution'):file_name.find('P')].replace('_',' ')
    p=float(file_name[file_name.find('P')+1:file_name.find('N')].replace('_','.'))
    N=int(file_name[file_name.find('N')+1:file_name.find('.pkl')])
    solutions[metodo,p,N]=dill.load(open(file_name,'rb') )
    ps.append(p)
    Ns.append(N)
    metodos.append(metodo)

ps=np.unique(ps)
Ns=np.unique(Ns)
metodos=list(set(metodos))
hs=[]
mallas={}
for N in Ns:
    mesh=Mesh('MallaN'+str(N)+'.xml')
    mallas[N]=mesh
    hs.append(mesh.hmax())
n_fig=0
n_cols=2
titutolos={'lenthg step':'Step Length','grad norm':r'$||DF||$','objective':r'$Objective function$',
          'lagrange multiplier':r'Lagrange multiplier $\mu$'}
for metodo in metodos:
    for i,p in enumerate(ps):
        n_fig+=1
        plt.figure(n_fig,figsize=(10,10))
        
        for j,N in enumerate(Ns):
            sol=solutions[metodo,p,N]
            n_subplots=len(titutolos)
            n_rows=int((n_subplots+1)/n_cols)
            for k,nombre in enumerate(titutolos):

                plt.subplot(n_rows,n_cols,k+1)
                plt.plot(sol[nombre],label='h={:.4}'.format(hs[j]))
                plt.title('Convergence history of '+titutolos[nombre])
                plt.ylabel(titutolos[nombre])
                plt.xlabel('iteration')
            
        for j in range(len(titutolos)):
            plt.subplot(n_rows,n_cols,j+1)
            plt.legend()
            plt.grid()
        plt.suptitle('Convergence history for p={:.2f}'.format(p)+' using '+metodo)
        plt.tight_layout()
        plt.savefig(save_dir+'ResumePLot'+metodo+str(p).replace('.','_')+'.png')
        plt.show()
for i,p in enumerate(ps):

    n_fig+=1
    plt.figure(n_fig,figsize=(5,5))
    for metodo in metodos:        
        objective=[ solutions[metodo,p,N]['objective'][-1] for N in Ns]
        
        plt.plot(hs,objective,'-*',label=metodo)
    plt.title('Convergence history of problem value for p='+str(p))
    plt.xlabel('Mesh size')
    plt.ylabel('Objective Function')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(save_dir+'MeshConvergenceRate'+str(p).replace('.','_')+'.png')
    
    plt.show()        

for p in ps:
    plt.figure(n_fig,figsize=(10,10))
    k=1
    n_fig+=1
    for i,metodo in enumerate(metodos):
        
        N=Ns[-1]
        sol=solutions[metodo,p,N]
        mesh=mallas[N]#Mesh('MallaN'+str(N)+'.xml')
        V,Theta= create_spaces(mesh)
        
        
        theta=Function(Theta)
        theta.vector()[:]=sol['theta vec']
        
        u=Function(V)
        u.vector()[:]=sol['u vec']
        
        plt.subplot(len(metodos),2,k)
        k+=1
        
        im=plot(u)
        plt.colorbar(im)
        plt.title('u obtained by '+metodo)
        plt.subplot(len(metodos),2,k)
        k+=1
        
        im=plot(theta)
        cbar=plt.colorbar(im)
        im.set_clim([0,1])
        plt.title(r'$\theta$ obtained by '+metodo)
        plt.xlabel('x')
        plt.ylabel('y')
    plt.suptitle('Solutions for p='+str(p))
    plt.tight_layout()
    plt.savefig(save_dir+'SolucionP'+str(p).replace('.','_')+'.png')

for p in ps:
    n_fig+=1
    plt.figure(n_fig,figsize=(7,5))
    tiempos={metodo : [solutions[metodo,p,N]['total time'] for N in Ns] for metodo in metodos}
    for metodo in tiempos:
        plt.scatter(hs,tiempos[metodo] ,label=metodo)
    plt.xlabel('Mesh size')
    plt.ylabel('CPU Time')
    plt.title('CPU Time per algorithm for p={:.2f}'.format(p))
    plt.legend()
    plt.grid(True)
    plt.savefig(save_dir+'Time'+str(p).replace('.','_')+'.png')
        
