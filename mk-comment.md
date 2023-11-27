Overall, I think the problem you were faced with was a numerical precision problem: 

In this function you had (rc/r) in the LJ potential, but it should be (1/r). The additional term (1/rc)^12 - (1/rc)^6 in the potential does not enter the force and can therefore be ignored. I introduced sigma=1 and used it in the force. The rc is the cutoff distance of the force, and not appearing anywhere else in the code, which is now the case: 
    
    def Lennard_Jones(r, r_c, k_LJ, norm_vect): #r = float, norm_vect from get_dist to get direction of force
        r, r_c = float(r), float(r_c)
        sigma = 1; # MK ADDED HERE AND BELOW 
        delta = 1e-5*sigma  # MK modified
    
        if r >= r_c or r == 0:
            return (0, 0)
    
        else:
            function_value = k_LJ * ((sigma/r)**12 - (sigma/r)**6) # MK modified
            gradient_value = (k_LJ * ((sigma/(r+delta))**12 - (sigma/(r+delta))**6) - function_value) / delta
            gradient = tuple(round(entry*gradient_value, 15) for entry in norm_vect)

        return gradient

Here you had r_c = 0.1, but this is too low, use r_c = 3.4 as written in the task description. 
        
      def get_LJ_forces(N, L, X, Y):
      k_LJ, r_c = 1, 3.4   # MK modified
      ...

Maybe there was a misinterpretation. If you use the values mentioned in the examples, like g = 3.5,
and T=0.3 instead of T=20 (T should not exceed T=2 or so) and k_S < 1 like k_S = 0.1 instead of k=10000, 
your code returns F=0 for the initial configuration. 

      def main():
      delete_folder_contents()
      # g, N, T, dt, k_S, dampening, steps = 0.1, 4, 20, 0.003, 1e4, 0.00, 1_000
      g, N, T, dt, k_S, dampening, steps = 3.5, 20, 0.3, 0.005, 0.1, 0.0, 1   # MK modified

Can you make a test using these parameters? Use parameter values mentioned in the task description and compare with the figures shown. And let me know immediately if you experience other problems. Overall, the code looks very good at first glance, but I didn't do any testing. 


      
