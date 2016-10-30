# inparam.nu
# created by Kuangdai on 28-Jun-2016 
# parameters needed to determine the integer field Nu(s, z)
# NOTE:
# a) angles are measured in degrees and distances in kilometers
# b) string-typed parameters (except file names) are case insensitive
# c) bool-typed parameters can be specified by 1/0, true/false, yes/no and on/off
# d) the prefix of input files is path_of_executable/input/, similar for output



# ================================== Type of Nu ==================================
# method to determine Nu field [constant / empirical / wisdom]
NU_TYPE                                     constant

# use fftw "lucky" numbers or not
# http://www.fftw.org/fftw2_doc/fftw_3.html
# FFTW is best at handling sizes of the form 2^a 3^b 5^c 7^d 11^e 13^f, 
# where e+f is either 0 or 1, and the other exponents are arbitrary.
# We call numbers of this form lucky numbers.
NU_FFTW_LUCKY_NUMBER                        true



# =================================== constant ===================================
# Nu(s, z) = NU_CONST
# NU_CONST = 2 is necessary and sufficient for 1D simulations
NU_CONST                                    2



# =================================== empirical ==================================
# use an empirical equation to determine Nu(s, z)
# Eqn. (69) in Kuangdai et al., Geophys. J. Int. (2016), Efficient global wave...  
# This empirical equation proves efficient and robust for global tomographic models.
# We suggest users simply increase/decrease NU_EMP_REF for better accuracy/performance. 
# A few trial simulations are aften needed for a certain combination of models. 

# basic reference value
NU_EMP_REF                                  24

# global minimum value
NU_EMP_MIN                                  8

# scale NU_EMP_REF by s-coordinate, i.e., distance to axis
# Fs = (s / R0) ^ NU_EMP_POW_AXIS
NU_EMP_SCALE_AXIS                           true
NU_EMP_POW_AXIS                             1.0

# scale NU_EMP_REF by epicentral distance (theta)
# F_theta = 1. + (NU_EMP_FACTOR_PI - 1.) * ((theta - NU_EMP_THETA_START) / (pi - NU_EMP_THETA_START)) ^ NU_EMP_POW_THETA
NU_EMP_SCALE_THETA                          true
NU_EMP_POW_THETA                            3.0
NU_EMP_FACTOR_PI                            5.0
NU_EMP_THETA_START                          45.0

# scale NU_EMP_REF by depth (km), to enhance surface wave resolution
# depth = 0 (surface)           => Fd = NU_EMP_FACTOR_SURF
# depth = NU_EMP_DEPHT_START    => Fd = NU_EMP_FACTOR_SURF
# depth = NU_EMP_DEPTH_END      => Fd = 1.0 
NU_EMP_SCALE_DEPTH                          true
NU_EMP_FACTOR_SURF                          2.0
NU_EMP_DEPTH_START                          200.0
NU_EMP_DEPTH_END                            300.0



# ================================== wisdom ==================================
# What is Wisdom? 
#     AxiSEM3D tries to learn the optimized Nu(s, z) in the next simulation,
#     and the result is written to a file, which can then be used repeatedly
#     in simulations with a similar scenario. Such a file is called a Wisdom,
#     a name borrowed form FFTW.
# How to make a Wisdom?
#     Learn a Wisdom by setting NU_WISDOM_LEARN true. The starting field, 
#     Nu_start(s, z), can be constant, empirical, or an existent Wisdom.   
#     The learnd field is always SMALLER than Nu_start(s, z). So, choose a 
#     Nu_start(s, z) larger than possibly required.  
# How to use a Wisdom?
#     Use a Wisdom by setting NU_TYPE to "wisdom", and specifying the Wisdom file
#     in NU_WISDOM_REUSE_INPUT. One may adjust the field by changing NU_WISDOM_REUSE_FACTOR.
#     Best practice: learn at a low frequency and reuse at higher frequencies.
#     You may NOT reuse a Wisdom if one of the following parameters changes a BIT:
#     a) 3D model, either volumetric or geometric
#     b) source location and depth
#     c) total record length

# If set true, AxiSEM3D will try to learn a Wisdom in the next simulation
# Wisdom learning slows the simulation down, but does not affect results 
NU_WISDOM_LEARN                             false

# balance accuracy/performance of the learned Wisdom [accurate / balance / fast]
NU_WISDOM_LEARN_AIM                         balance   

# interval of Wisdom learning
NU_WISDOM_LEARN_INTERVAL                    100 

# file to save the learned Wisdom
NU_WISDOM_LEARN_OUTPUT                      name.nu_wisdom

# Wisdom that will be used in the next simulation 
NU_WISDOM_REUSE_INPUT                       name.nu_wisdom

# factor multiplied to Nu(s, z) specified in NU_WISDOM_REUSE_INPUT
NU_WISDOM_REUSE_FACTOR                      1.0
