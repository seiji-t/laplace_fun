import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    args = sys.argv

    _, nmax, nt = args

    azi_values = np.linspace(0,360, nt)
    for it in range(nt):
        azi = azi_values[it]
        fig = plt.figure(figsize=(12,12))
        kount=1
        for n in range(nmax):
            for m in range(n+1):
                data = np.loadtxt('scripts/data/{:02d}_{:02d}_{:04d}.txt'.format(n,m,it), delimiter=',')
                n_mu = 64
                n_theta = 128

                theta = data[:,0].reshape(n_mu, n_theta)
                mu = data[:,1].reshape(n_mu, n_theta)
                u = data[:,2].reshape(n_mu, n_theta)

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

def nn(n):
    if n == 0:
        return 1
    return n * nn(n-1)
def nn(l,m):
    return math.sqrt( (2*l+1) * math.factorial(l-m)/2/math.pi/math.factorial(l+m) )

nn(0,0)