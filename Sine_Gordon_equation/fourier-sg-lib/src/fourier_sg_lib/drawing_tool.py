import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
import numpy as np

def generate_gif(x,U_t,real_time,filename="soliton.gif",
                 interval=50,U_a=[],title='', y_min=-7, y_max =  7,
                 msr = False,Umsr=[],xmsr=[],
                 mark=False,Umark=[]):
    if filename[-4:]!=".gif":
        filename += ".gif"

    fig, ax = plt.subplots()
    line, = ax.plot(x, U_t[0],label='numerical soultion')
    if msr:
        ax.scatter(xmsr,Umsr,label='measuring points',s=6,color='red')
    if len(U_a)!=0:
        line_a, = ax.plot(x,U_a[0],label='analytic solution',linestyle = 'dotted')
    if mark:
        line_mark = ax.axvline(x=Umark[0],ymin=-5,ymax=5,color='red')
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("x")
    ax.set_ylabel("U(x)")
    fig.suptitle(title)
    ax.legend()

    def update(frame):
        line.set_ydata(U_t[frame])
        if len(U_a)>0:
            line_a.set_ydata(U_a[frame])
        if mark:
            line_mark.set_xdata([Umark[frame],Umark[frame]])
        ax.set_title(f"Time {real_time[frame]:.2f}")
        return line,

    anim = FuncAnimation(fig, update, frames=U_t.shape[0]-1, interval=interval, blit=True)


    anim.save(filename, writer=PillowWriter(fps=1000//interval))
    plt.close(fig)

    print(f"GIF saved as '{filename}'")


def giff_comp(x,U_t,real_time,labels,title,filename="soliton.gif",interval=50):
    if filename[-4:]!=".gif":
        filename += ".gif"

    y_min, y_max = -7, 7

    fig, ax = plt.subplots()
    lines = []
    for U,label in zip (U_t,labels):
        line, = ax.plot(x, U[0],label=label)
        lines.append(line)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("x")
    ax.set_ylabel("U(x)")
    fig.suptitle(title)
    ax.legend(loc='upper right')

    def update(frame):
        for U, line in zip(U_t,lines): 
            line.set_ydata(U[frame])
        ax.set_title(f"Time {real_time[frame]:.2f}")
        return line,

    anim = FuncAnimation(fig, update, frames=U_t[0].shape[0]-1, interval=interval, blit=True)


    anim.save(filename, writer=PillowWriter(fps=1000//interval))
    plt.close(fig)

    print(f"GIF saved as '{filename}'")

def plot_a_frame(x,U_t,real_time,nframe,labels,title,filename,ls,ylim=[-7,7]):
    if filename[-4:]!=".png":
        filename += ".png"
    y_min, y_max = ylim[0],ylim[1] 
    
    fig, ax = plt.subplots()
    lines = []
    for U,label,l in zip (U_t,labels,ls):
        line, = ax.plot(x, U[nframe],label=label,linestyle=l)
        lines.append(line)
    ax.set_title(f"Time {real_time[nframe]:.2f}")
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("x")
    ax.set_ylabel("U(x)")
    fig.suptitle(title)
    ax.legend()
    fig.savefig(filename)
    print(f"Plot {filename} saved!")


def plot_3d(x,U_t,real_time,filename="soliton.png",el=20,azm=15):
    if filename[-4:]!=".png":
        filename += ".png"
    fig= plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(projection='3d')
    lines = []

    for U in U_t:
        X, Tmesh = np.meshgrid(x, real_time)
        ax.plot_surface(Tmesh,X,U,rstride=2, cstride=2, cmap='viridis')
        ax.view_init(elev=el, azim=azm) 
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.set_zlabel("U(x)")
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
    ax.set_position([0, 0, 1, 1])
    ax.set_box_aspect([1,1.33,1])
    fig.savefig(filename)
    print(f"Plot {filename} saved!")

def plot_cm(x,U_t,real_time,filename="soliton.png",two_breahters=False,U_marker=[],xlim=[]):
    
    fig, ax = plt.subplots()
    if len(xlim)>0:
        ax.set_xlim(xlim)
    for U in U_t:
        pcm = ax.pcolormesh(
                x,
                real_time,
                U,
                cmap='viridis',
                shading='auto'
                )
        fig.colorbar(pcm, ax=ax, label='U(x,t)')
    if two_breahters:
        
        ax.plot(U_marker,real_time,linestyle='--',color="gray",label = "Trajectory of uninteracting breather")
    ax.set_xlabel('x')
    ax.set_ylabel('time')
    ax.legend()
    fig.savefig(filename)
    print(f"Plot {filename} saved!")

def plt_eps(epsm, epsi,y_mass,y_inter,title,y_label):
    fig, ax = plt.subplots()
        
    ax.plot(epsm,y_mass,c='lightgrey',linestyle='--')
    ax.plot(epsi,y_inter,c='grey',linestyle='--')
        
    ax.scatter(epsm,y_mass,s=7,label= f"{y_label} mass perturbation")
    ax.scatter(epsi,y_inter,s=7,marker='^',label= f"{y_label} interactions perturbation")
        
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel("$\epsilon$")
    
    ax.legend()
    ax.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))
    ax.set_xscale('symlog', linthresh=1e-7)
    filename = f"eps_{y_label}.png"
    fig.savefig(filename)
    print(f"Plot {filename} saved!")