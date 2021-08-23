import log_gaussian_boost_python 
import time
import numpy as np 

np.random.seed(1234567890)

#-8445778 +/- 2822
#-8448073.607605366

#(lst-dev_l) > python log_gaussian.py 
#tottime/nPoints  =  0.11108147382736205
#(lst-dev_l) > python log_gaussian_numba_jit.py 
#jit (first call)  = 0.5462789535522461
#jit (second call) = 0.12479305267333984
#tottime/nPoints   = 0.11998658180236817
#(lst-dev_l) > python log_gaussian_numba_CC.py 
#CC (first call)  = 0.09273552894592285
#CC (second call) = 0.09696602821350098
#tottime/nPoints  =  0.09173028707504273
#(boost_python) > python log_gaussian_boost_python.py 
#tottime/nPoints  =  0.10884714603424073

#python    0.11    100%
#numba jit 0.12    110%
#numba CC  0.09     82%
#boost     0.11    100%

def main():

    #x_shape     =  (196, 609, 1)
    #mean_shape  =  (196, 609, 38)
    #sigma_shape =  (196, 609, 38)
    #i_n = 2
    #j_n = 3
    #k_n = 4
    i_n = 196
    j_n = 609
    k_n = 38
    x_shape     =  (i_n, j_n, 1)
    mean_shape  =  (i_n, j_n, k_n)
    sigma_shape =  (i_n, j_n, k_n)
    #
    x = np.random.normal(3, 1.0, size=x_shape)
    mean = np.random.normal(3, 1.0, size=mean_shape)
    sigma = np.random.normal(2.0, 0.001, size=sigma_shape)
    log_pdf = np.zeros(shape=sigma_shape)

    #print("i"," ","j"," ","k"," ","val")
    #for i in range(i_n) :
    #    for j in range(j_n) :
    #        for k in range(k_n) :
    #            print(i," ",j," ",k," ",mean[i][j][k])
    
    nPoints = 100
    tottime = 0
    
    for i in range(nPoints):
        time01 = time.time()
        log_gaussian_boost_python.log_gaussian(x,mean,sigma,log_pdf)
        time02 = time.time()
        delta=(time02-time01)
        tottime = tottime + delta

    print('tottime/nPoints  = ',tottime/nPoints)


if __name__ == "__main__":
    main()
