# Soils Animated: Part 1

Animation is a powerful tool for demonstrating concepts. It is engaging, makes abstract scientific ideas easier to understand, and effectively communicates ideas to a diverse audience.

In this blog post, I will demonstrate a simple soil dynamics model (Laio et al., 2001; Rodriguez-Iturbe et al., 1999) in animation. I like this model and have been using in my own research, because it helps us conceptualize soil hydrologic system at the field scale with a small number of key variables. This 'toy model' for ecohydrologists have been used for many interesting experiments: D'Odorico and Porporato (2004) explained soil moisture seasonality, and Entekhabi and Rodriguez-Iturbe (1994) investigated the impacts of spatio-temporal aggregation on characterizing soil dynamics heterogeneity, and many others. In this blog, I will try to illustrate the basic concepts of the model, with Part 1 focusing on how soil moisture dynamics can be translated in different spaces using this model. 

## Soil Moisture Loss Function to Describe the Dynamics
The model is called the soil moisture loss function, that relates the rate of loss from soil volume (-dtheta/dt) at a given soil moisture level (theta). It is a type of Storage-Flux relation model in Hydrology, describing how Fluxes (in this case, drainage and evapotranspiration) are regulated by Storage (soil moisture) availability. This function describes soil moisture dynamics in two major stages: Drainage and Evapotranspiration (for simplicity, hygroscopic water are neglected here).

### Drainage Stage 
Think about a rainy day. Shortly after rainfall ceases, soil becomes saturated or extremely wet. If you step on it, soil gets muddy and squishy. In this state, soil macropores are filled with water which can drain out quickly due to gravity. Some water might run off the soil surface. The first stage of the loss function describes this state, and the processes are expressed as an exponential function. 

![Alt text](.\out\drainage.gif "Drainage")


### Evapotranspiration Stage
After a day or few days, the soil becomes semi-dry - your back might still get wet if you sit on it, but when touched, the soil is moist without water immediately seeping out. In this state, water in the soil macropores are mostly drainagd out. Much of the water is held against gravity within the soil micropores because the soil suction force.

Flux coming out of this semi-dry soil is dominated by evapotranspiration (ET); water are be pulled out from soil to the atmosphere mainly through the ET process. When the soil is wet, the transpiration process reaches its maximum rate with fully opened plants stomata, and ET can reach its maximum potential (note: though evaporation can continue independently; Krell et al., 2021). As soil dries out, a threshold called critical soil moisture theta* is reached, where plants begin to experience water stress. Below that point, ET decreases proportionally with soil moisture level (although this can be nonlinear; Araki et al., 2024). This ET dynamics is represented by the piecewise-linear function in the loss function space. 

![Alt text](.\out\ET.gif "ET")

## Usefulness of the Soil Moisture Loss Function
This loss function model treates soil volume as a system responding to pulse inputs of rainfall (which I plan to elaborate on and animate in Part 2). The soil responds in a prescribed manner based on the combination of soil physical properties, external forcing, and current system state. The loss function describes this underlying system as shown in the two animations shown earlier: the left panel, the loss function, describes the underlying system. 

Because the loss function is an ordinary differential equation (ODE), we can move our perspectives between the underlying system space and the observable space. ODE can be numerically or analytically solved, and the analytical solution of the loss function (left panel) is the drydown curve (right panel). The soil moisture drying behavior is something we can observe using soil moisture sensors. Therefore, from the observed drydown curve, we can back-estimate the underlying system, the loss function shape. 

One of the advantages of being able to move between different spaces is that it reduces uncertainties in model estimation. Here's why: 

The loss function space describes more general patterns of the system but is prone to data error when applied to observation data. Both y-variable and x-variable contain theta which is subject to observation uncertainty. In the figure below, you can see how loss function estimates are affected observation errors in the drydown curves (in the animation, gaussian noises are artifically generated). This situation calls for bivariate analysis but they are more difficult to implement. 

------------------- Add animation of data error mapped with loss function model ---------------------------------

However, in the drydown space (which is the analytical solution of the loss function), the x-variable becomes time - a quantity without observation uncertainty. In this case, only the y-variable has uncertainty, allowing the use of more straightforward nonlinear least squares fitting. Therefore, although the observed drydown may not cover capture the entire range of the loss function, the drydown models are advantageous in mapping the dynamics, as they allow us to deal observation errors only in the x-variable (which I assume is more robust). 

------------------- Add animation of data error mapped with drydown model ---------------------------------

## Power of the Animation
Science requires abstraction of processes, and we use models to capture these processes. With the use of animation, we can easily visualize the model behavior, and sometimes visualize the more abstracted spaces.  (shameless plug),Animation helps us to wrap around our head to move between observation space and theoretical spaces, just like the drydown space and the loss function space. 

What's great about this DeepGroundwater platform is that it allows us to share animations (a shameless plug)! While animations are difficult to include in traditional scientific formats like PDFs, this platform freely accommodates images, animations, videos, and even JavaScript animations to share ideas. We welcome scientific & artistic contributions here.

Part 2 of this blog post will demonstrate the challenges of using this soil moisture loss function model - which is actullly a bucket type of model that was pointed out as an oversimplification of the reality in [our first blogpost](https://deepgroundwater.com/blog/hydrology-is-flat-and-its-buckets-all-the-way-down/). I plan to describe how loss function plays out in longer soil moisture timeseries with multiple pulse rainfall input, and how we may or maynot be able to derive information from it (with animations, of course!). Stay tuned.

## Code Availability 
I created these animations using Python's matplotlib, which you can find here: https://github.com/RY4GIT/drydown-viz.

## Reference & Acknolwegements
The ideas of this blog post comes from a lots of different pieces of ideas discussed with my advisors and collaborators, thank you all! 



Laio, F., Porporato, A., Fernandez-Illescas, C. P., & Rodriguez-Iturbe, I. (2001). Plants in water-controlled ecosystems: Active role in hydrologic processes and responce to water stress IV. Discussion of real cases. Advances in Water Resources, 24(7), 745–762. https://doi.org/10.1016/S0309-1708(01)00007-0

Rodriguez-Iturbe, I., Porporato, A., Ridolfi, L., Isham, V., & Coxi, D. R. (1999). Probabilistic modelling of water balance at a point: the role of climate, soil and vegetation. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences. https://doi.org/10.1098/rspa.1999.0477

D’Odorico, P., & Porporato, A. (2004). Preferential states in soil moisture and climate dynamics. Proceedings of the National Academy of Sciences of the United States of America, 101(24), 8848–8851. https://doi.org/10.1073/pnas.0401428101

Entekhabi, D., & Rodriguez-Iturbe, I. (1994). Analytical framework for the characterization of the space-time variability of soil moisture. Advances in Water Resources, 17(1), 35–45. https://doi.org/10.1016/0309-1708(94)90022-1

Krell, N. T., Morgan, B. E., Gower, D., & Caylor, K. K. (2021). Consequences of dryland maize planting decisions under increased seasonal rainfall variability. Water Resources Research, 57(9). https://doi.org/10.1029/2020wr029362

Araki, R., Morgan, B., McMillan, H. K., & Caylor, K. (2024). Nonlinear soil moisture loss function reveals vegetation responses to water availability. Authorea Preprints. https://doi.org/10.22541/essoar.172251989.99347091
