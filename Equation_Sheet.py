import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed
import Equation_Sheet as es
gamma=4.5e-8

def derivs2(ivec,t,M,S):
        x=ivec[0]
        vx=ivec[1]
        y=ivec[2]
        vy=ivec[3]                                                           #this is the star x,y positions and velocities
        
        X=ivec[4]
        Vx=ivec[5]
        Y=ivec[6]
        Vy=ivec[7]                                                           #this is the galaxy x,y positions and velocities
        
        dx=vx
        dy=vy
        
        dX=Vx
        dY=Vy
        
        r=np.sqrt(x**2+y**2)
        R=np.sqrt(X**2+Y**2)      
        
        rho=np.sqrt((X-x)**2+(Y-y)**2)
        rhoy=Y-y
        rhox=X-x
        
        dvx = -gamma*((M/r**3)*x-(S/rho**3)*rhox+(S/R**3)*X)
        dVx = -gamma*((M+S)/R**3)*X
        
        dvy = -gamma*((M/r**3)*y-(S/rho**3)*rhoy+(S/R**3)*Y)
        dVy = -gamma*((M+S)/R**3)*Y                                          #partial differential equations
        return np.array([dx,dvx,dy,dvy,dX,dVx,dY,dVy])
    
def sol(ic, max_t, ntimes, M, S):
    t = np.linspace(0.0, max_t, ntimes)
    s = odeint(derivs2, ic, t, args=(M,S))                                   #solves the partial differential equations
    return s

def generate_star_ics(M, radii, nstars): 
    """generates initial conditions for n number of stars at any                                                                               radii for galaxy M (radii and nstars must be passed as lists 
    or arrays of the same length)"""
    ics = []
    for i in range(len(radii)):
        r = radii[i]
        N = nstars[i]
        theta = np.arange(0.0, 2.0*np.pi, 2.0*np.pi/N)
        for t in range(N):
            x = r*np.cos(theta[t])
            y = r*np.sin(theta[t])
            v = np.sqrt((gamma*M)/r)
            vx = -v*np.sin(theta[t])
            vy = v*np.cos(theta[t])
            ics.append(np.array([x,vx,y,vy]))
    return ics

def angle(Y,r0):
    return np.arctan(2.0*r0/Y) #generates angles for generate_gala_ics

def generate_gala_ics(M,S,r0,m):
    """Forms the moving galaxies' (S) initial conditions"""
    Y = m*r0
    X = (-Y**2.0)/(4.0*r0)+r0
    V=np.sqrt(2.0*gamma*(M+S)/np.sqrt(X**2+Y**2))
    Vx = V*np.cos(angle(Y,r0))
    Vy = -V*np.sin(angle(Y,r0))
    return np.array([X,Vx,Y,Vy])

def make_stars(M,S,radii,nstars,r0,m,max_t,ntimes): 
    """takes the star and galaxy initial conditions and returns them for the ploting function to reduce the lag on interact"""
    star_ics = generate_star_ics(M,radii,nstars) # [np.array([x,vx,y,vy]), np.array([]), ...]
    gala_ics = generate_gala_ics(M,S,r0,m) # np.array([X,Vx,Y,Vy])
    #ics[0] -> x, vx, y, vy
    starx = []
    stary = []
    for sic in star_ics:
        ic = np.hstack([sic, gala_ics])
        result = sol( ic, max_t, ntimes, M, S)
        starx.append(result[:,0])
        stary.append(result[:,2])
        galax = result[:,4]
        galay = result[:,6]
    starx = np.transpose(np.array(starx))
    stary = np.transpose(np.array(stary))
    return starx, stary, galax, galay

def plot_solution(starx, stary, galax, galay, j, lim):
    """plots the solution recived from the make_stars equation"""
    px=np.linspace(-100,100,100)
    r0=25.0
    py=-px**2/(4.0*r0)+r0
    plt.plot(py,px,color='orchid')
    plt.scatter(starx[j],stary[j],color='b')
    plt.scatter(galax[j],galay[j],color='lime')
    plt.scatter(0,0,color='r')
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    
def tiles():
    """makes the 4x4 subplots for each question"""
    starx, stary, galax, galay=es.make_stars(1e10,1e10,[4,8,12,16,20],[13,19,24,30,36],25,3,50,500)
    plt.figure(figsize=(25,25))
    plt.subplot(4,4,1)
    es.plot_solution(starx,stary,galax,galay,100,40)
    plt.subplot(4,4,2)
    es.plot_solution(starx,stary,galax,galay,135,40)
    plt.subplot(4,4,3)
    es.plot_solution(starx,stary,galax,galay,145,40)
    plt.subplot(4,4,4)
    es.plot_solution(starx,stary,galax,galay,150,40)
    plt.subplot(4,4,5)
    es.plot_solution(starx,stary,galax,galay,170,40)
    plt.subplot(4,4,6)
    es.plot_solution(starx,stary,galax,galay,200,40)
    plt.subplot(4,4,7)
    es.plot_solution(starx,stary,galax,galay,230,40)
    plt.subplot(4,4,8)
    es.plot_solution(starx,stary,galax,galay,250,40)
    plt.subplot(4,4,9)
    es.plot_solution(starx,stary,galax,galay,275,40)
    plt.subplot(4,4,10)
    es.plot_solution(starx,stary,galax,galay,300,40)
    plt.subplot(4,4,11)
    es.plot_solution(starx,stary,galax,galay,330,40)
    plt.subplot(4,4,12)
    es.plot_solution(starx,stary,galax,galay,350,40)
    plt.subplot(4,4,13)
    es.plot_solution(starx,stary,galax,galay,370,40)
    plt.subplot(4,4,14)
    es.plot_solution(starx,stary,galax,galay,400,40)
    plt.subplot(4,4,15)
    es.plot_solution(starx,stary,galax,galay,450,40)
    plt.subplot(4,4,16)
    es.plot_solution(starx,stary,galax,galay,499,40)

def generate_star_ics_dir(M, radii, nstars):
    """same as avove make_stars function but swiches the direction of the orbit"""
    ics = []
    for i in range(len(radii)):
        r = radii[i]
        N = nstars[i]
        theta = np.arange(0.0, 2.0*np.pi, 2.0*np.pi/N)
        for t in range(N):
            x = r*np.cos(theta[t])
            y = r*np.sin(theta[t])
            v = np.sqrt((gamma*M)/r)
            vx = v*np.sin(theta[t])
            vy = -v*np.cos(theta[t])
            ics.append(np.array([x,vx,y,vy]))
    return ics

def make_stars_dir(M,S,radii,nstars,r0,m,max_t,ntimes):
    star_ics = generate_star_ics_dir(M,radii,nstars)
    gala_ics = generate_gala_ics(M,S,r0,m) # np.array([X,Vx,Y,Vy])
    #ics[0] -> x, vx, y, vy
    starx = []
    stary = []
    for sic in star_ics:
        ic = np.hstack([sic, gala_ics])
        result = sol( ic, max_t, ntimes, M, S)
        starx.append(result[:,0])
        stary.append(result[:,2])
        galax = result[:,4]
        galay = result[:,6]
    starx = np.transpose(np.array(starx))
    stary = np.transpose(np.array(stary))
    return starx, stary, galax, galay
    """generates initial conditions for n number of stars at any                                                                               radii for galaxy M (radii and nstars must be passed as lists 
    or arrays of the same length) for the direct passage"""

def tiles_dir():
    """same as above tiles function but for a direct passage"""
    starx, stary, galax, galay=es.make_stars_dir(1e10,1e10,[4,8,12,16,20],[13,19,24,30,36],25,3,50,500)
    plt.figure(figsize=(25,25))
    plt.subplot(4,4,1)
    es.plot_solution(starx,stary,galax,galay,100,40)
    plt.subplot(4,4,2)
    es.plot_solution(starx,stary,galax,galay,135,40)
    plt.subplot(4,4,3)
    es.plot_solution(starx,stary,galax,galay,145,40)
    plt.subplot(4,4,4)
    es.plot_solution(starx,stary,galax,galay,150,40)
    plt.subplot(4,4,5)
    es.plot_solution(starx,stary,galax,galay,170,40)
    plt.subplot(4,4,6)
    es.plot_solution(starx,stary,galax,galay,200,40)
    plt.subplot(4,4,7)
    es.plot_solution(starx,stary,galax,galay,230,40)
    plt.subplot(4,4,8)
    es.plot_solution(starx,stary,galax,galay,250,40)
    plt.subplot(4,4,9)
    es.plot_solution(starx,stary,galax,galay,275,40)
    plt.subplot(4,4,10)
    es.plot_solution(starx,stary,galax,galay,300,40)
    plt.subplot(4,4,11)
    es.plot_solution(starx,stary,galax,galay,330,40)
    plt.subplot(4,4,12)
    es.plot_solution(starx,stary,galax,galay,350,40)
    plt.subplot(4,4,13)
    es.plot_solution(starx,stary,galax,galay,370,40)
    plt.subplot(4,4,14)
    es.plot_solution(starx,stary,galax,galay,400,40)
    plt.subplot(4,4,15)
    es.plot_solution(starx,stary,galax,galay,450,40)
    plt.subplot(4,4,16)
    es.plot_solution(starx,stary,galax,galay,499,40)

def make_stars_heavy(M,S,radii,nstars,r0,m,max_t,ntimes):
    star_ics = generate_star_ics(M,radii,nstars) # [np.array([x,vx,y,vy]), np.array([]), ...]
    gala_ics = generate_gala_ics(M,S,r0,m) # np.array([X,Vx,Y,Vy])
    #ics[0] -> x, vx, y, vy
    starx = []
    stary = []
    for sic in star_ics:
        ic = np.hstack([sic, gala_ics])
        result = sol( ic, max_t, ntimes, M, S)
        starx.append(result[:,0])
        stary.append(result[:,2])
        galax = result[:,4]
        galay = result[:,6]
    starx = np.transpose(np.array(starx))
    stary = np.transpose(np.array(stary))
    return starx, stary, galax, galay

def tiles_heavy():
    """same as above tiles function but for a heavy mass passage"""
    starx, stary, galax, galay=es.make_stars_heavy(1e10,1e15,[4,8,12,16,20],[13,19,24,30,36],25,3,.25,500)
    plt.figure(figsize=(25,25))
    plt.subplot(4,4,1)
    es.plot_solution(starx,stary,galax,galay,100,40)
    plt.subplot(4,4,2)
    es.plot_solution(starx,stary,galax,galay,135,40)
    plt.subplot(4,4,3)
    es.plot_solution(starx,stary,galax,galay,145,40)
    plt.subplot(4,4,4)
    es.plot_solution(starx,stary,galax,galay,150,40)
    plt.subplot(4,4,5)
    es.plot_solution(starx,stary,galax,galay,170,40)
    plt.subplot(4,4,6)
    es.plot_solution(starx,stary,galax,galay,200,40)
    plt.subplot(4,4,7)
    es.plot_solution(starx,stary,galax,galay,230,40)
    plt.subplot(4,4,8)
    es.plot_solution(starx,stary,galax,galay,250,40)
    plt.subplot(4,4,9)
    es.plot_solution(starx,stary,galax,galay,275,40)
    plt.subplot(4,4,10)
    es.plot_solution(starx,stary,galax,galay,300,40)
    plt.subplot(4,4,11)
    es.plot_solution(starx,stary,galax,galay,330,40)
    plt.subplot(4,4,12)
    es.plot_solution(starx,stary,galax,galay,350,40)
    plt.subplot(4,4,13)
    es.plot_solution(starx,stary,galax,galay,370,40)
    plt.subplot(4,4,14)
    es.plot_solution(starx,stary,galax,galay,400,40)
    plt.subplot(4,4,15)
    es.plot_solution(starx,stary,galax,galay,450,40)
    plt.subplot(4,4,16)
    es.plot_solution(starx,stary,galax,galay,499,40)
    
def tiles_lite():
    """same as above tiles function but for a light mass passage"""
    starx, stary, galax, galay=es.make_stars(1e10,5e9,[4,8,12,16,20],[13,19,24,30,36],25,3,50,500)
    plt.figure(figsize=(25,25))
    plt.subplot(4,4,1)
    es.plot_solution(starx,stary,galax,galay,100,40)
    plt.subplot(4,4,2)
    es.plot_solution(starx,stary,galax,galay,135,40)
    plt.subplot(4,4,3)
    es.plot_solution(starx,stary,galax,galay,145,40)
    plt.subplot(4,4,4)
    es.plot_solution(starx,stary,galax,galay,150,40)
    plt.subplot(4,4,5)
    es.plot_solution(starx,stary,galax,galay,170,40)
    plt.subplot(4,4,6)
    es.plot_solution(starx,stary,galax,galay,200,40)
    plt.subplot(4,4,7)
    es.plot_solution(starx,stary,galax,galay,230,40)
    plt.subplot(4,4,8)
    es.plot_solution(starx,stary,galax,galay,250,40)
    plt.subplot(4,4,9)
    es.plot_solution(starx,stary,galax,galay,275,40)
    plt.subplot(4,4,10)
    es.plot_solution(starx,stary,galax,galay,300,40)
    plt.subplot(4,4,11)
    es.plot_solution(starx,stary,galax,galay,330,40)
    plt.subplot(4,4,12)
    es.plot_solution(starx,stary,galax,galay,350,40)
    plt.subplot(4,4,13)
    es.plot_solution(starx,stary,galax,galay,370,40)
    plt.subplot(4,4,14)
    es.plot_solution(starx,stary,galax,galay,400,40)
    plt.subplot(4,4,15)
    es.plot_solution(starx,stary,galax,galay,450,40)
    plt.subplot(4,4,16)
    es.plot_solution(starx,stary,galax,galay,499,40)