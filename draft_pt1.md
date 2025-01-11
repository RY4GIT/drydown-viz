# Soil Volume as a Bucket: Part 1

Animation is a powerful tool for demonstrating concepts. In this blog post, I will demonstrate a simple soil dynamics model (Laio et al., 2001; Rodriguez-Iturbe xxxx) that I have been using in my own research, and illustrate the movement between different conceptual spaces using animations.

## Soil moisture loss function
The concept I want to describe here is the simplified soil moisture dynamics, called the soil moisture loss function, that relates the loss from soil volume at a given soil moisture level. It is a type of Storage-Flux relationship in Hydrology, describing how Fluxes (in this case, drainage and evapotranspiration) are regulated by Storage (soil moisture) levels. This function primarily describes soil moisture dynamics in two main stages : Drainage and Evapotranspiration (for simplicity, the dynamics after hygroscopic water are neglected here).

### Drainage stage 
After rainfall, soil becomes saturated or extremely wet. The soil gets muddy and squishy. Water is held in the macropores of the soil and can drain out quickly due to gravity. Some water might run off the surface.The first stage of the loss function describes this state and processes as an exponential (checkxxxx) function. 

![Alt text](.\out\drainage.gif "Drainage")


### Evapotranspiration stage
After a day or few days, the soil becomes semi-dry - your back might still get wet if you sit on it, but when touched, the soil is moist without water immediately seeping out. The water is held against gravity within the soil because the sorptive suction of soil micropores retains the water.

Flux coming out of this semi-dry soil is dominated by evapotranspiration (ET). This water can be pulled out from soil mainly through the ET process. When the soil is wet, the transpiration process reaches its maximum rate with fully opened plants stomata, and ET can reach its maximum potential (note: though evaporation can continue independently (Krell et al., 2021)). As soil dries out, a threshold called critical soil moisture is reached, where plants begin to experience water stress. Below that point, ET decreases proportionally with soil moisture level (although this can be nonlinear; Araki xxxxx, 2005). This ET dynamics is represented by the piecewise-linear function in the loss function space. 

![Alt text](.\out\ET.gif "ET")

# Usefulness of loss function
This loss function model treates soil volume as a system responding to pulse inputs of rainfall (which I plan to animate in Part 2). The soil responds in a prescribed manner based on its physical properties, external forcing, and current system state. The loss function is useful because it describes this underlying system that governs how soil volume responds to input pulses, as shown in the animation: the left panel, the loss function, describes the general system, and the right panel, describes the observation of soil mositure reponse to pulse rainfall.

Furthermore, the loss function is an ordinary differential equation (ODE) that allows us to move between drydown space and loss function space. One of its advantages lies in reducing uncertainties in model estimation. Here's why: 

The loss function space describes more general patterns of the system but is prone to data error when applied to observation data. Both y-variable and x-variable contain theta (soil moisture), which is subject to data uncertainty. In the figure below, you can see how loss function estimates are affected by fluctuations mimicing observation errors in the drydown curves. This situation calls for bivariate analysis but they are more difficult to implement (is it true?). 

------------------- Add animation of data error mapped with loss function model ---------------------------------

However, in the drydown space (i.e., which is the analytical solution of the loss function), the x-variable becomes time - a quantity without observation uncertainty. In this case, only the y-variable has uncertainty, allowing the use of more general nonlinear least squares analysis. Therefore, although the observed drydown may not cover capture the entire range of the loss function, the drydown models are advantageous in mapping the dynamics, as they can treat errors as Gaussian noise in the optimization algorithms. The estiamted optimal drydown models are then can be back-mapped to the more general, loss function space. 

------------------- Add animation of data error mapped with drydown model ---------------------------------

# Power of animation
Science requires abstraction of processes, and we use models to reproduce these processes. With the use of animation, we can easily moving between these spaces, and sometimes between observation space and theoretical spaces (just like the drydown space and the loss function space). I created these animations using Python's matplotlib, which you can find here: https://github.com/RY4GIT/drydown-viz.

Additionally (shameless plug), what's great about this DeepGroundwater platform is that it allows us to share animations! While animations are difficult to include in traditional scientific formats like PDFs, this platform freely accommodates images, animations, videos, and even JavaScript animations to share ideas. We welcome scientific & artistic contributions here.

Part 2 of this blog post will demonstrate the challenges of using the loss function to describe and derive information from observed soil drying behaviors (with animations, of course!). Stay tuned.

# Reference & Acknolwegements
The ideas of this blog post comes from a lots of different pieces of ideas discussed with my advisors and collaborators, thanks all! 


<!-- # Chalalnges of loss function -->
<!-- However, this is much of a bucket model that simplifies very hetegronegous conditions (see the blog post). In my experience, this model works pretty ok for point-based data. You can actually see the point of inflection from drainage stage to ET stage in observed drydown curves (this animation, again for the Part 2) and make sense of the parameters, as shown in previous studies. I hypothesize that this bucket type model may work for a small soil volume and provides an emergent pattersn of soil dynamics: emergent patterns that is aggregated from the core-scale equation, like Richards equation. However, there is also an influence of deepeer soil layer, that derives the pressure gradient, that actually the drives soil moisture dynamics. I am reading right now this book, also analysing the data, and making animations to make sense of all the dynamics. However, as shown in the curve on the right panel, the transition between different ET phases in the loss function is very difficult to distinguish in the drydown space. In fact, most studies determining critical soil moisture content do not derive this from the drydown space. Instead, they typically use either explore this parameter in the loss function space (which can be vulnerable to data uncertainty, as discussed in the next section) or combine soil moisture observations with ET observation from flux towers (xxxx). -->
