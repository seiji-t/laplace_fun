import matplotlib.pyplot as plt
import numpy as np
import sys
import csv

def read_eigenfunction(n,m):
    filename = 'data/{:02d}_{:02d}.txt'.format(n,m)
    theta = []
    mu = []
    u = []
    params = {}
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file, delimiter=',')
        for n,row in enumerate(csv_reader):
            if n == 0:
                params['omega'] = float(row[0])
                params['c'] = float(row[1])
                params['radius'] = float(row[2])
                params['n_mu'] = int(row[3])
                params['n_theta'] = int(row[4])
            else:
                theta.append(float(row[0]))
                mu.append(float(row[1]))
                u.append(float(row[2]))
    return params, theta, mu, u

def main():
    args = sys.argv

    _, nmax, nt = args
    nt = int(nt)
    nmax = int(nmax)
    
    azi_values = np.linspace(0,360, nt) * 0
    t_values = np.linspace(0, 2*np.pi, nt)

    for it, (t, azi) in enumerate(zip(t_values, azi_values)):
        print(t)
        fig = plt.figure(figsize=(12,12))
        kount=1
        for n in range(nmax):
            for m in range(n+1):
                params, theta, mu, u = read_eigenfunction(n,m)
                n_mu = params['n_mu']
                n_theta = params['n_theta']
                omega = params['omega']
                c = params['c']
                r = params['radius']

                theta = np.array(theta).reshape(n_mu, n_theta)
                mu = np.array(mu).reshape(n_mu, n_theta)
                u = (np.array(u).reshape(n_mu, n_theta) - r)*np.cos(omega * t/c) + r

                phi = np.arccos(mu)

                r = u
                x = r * np.sin(phi) * np.cos(theta)
                y = r * np.sin(phi) * np.sin(theta)
                z = r * np.cos(phi)
                ax = fig.add_subplot(nmax,nmax,kount, projection='3d')
                ax.plot_wireframe(x, y, z, color='black', linewidth=0.5)

                ax.set_box_aspect([1,1,1])
                ax.axis('off')
                ax.view_init(elev=30, azim=azi)
                plt.title(f"m={m}, n={n}", fontsize=16)
                kount+=1
            kount+= nmax-m-1

        plt.savefig("fig/fig_{:04d}.png".format(it))
        plt.close()    
        
if __name__ == "__main__":
    main()