Run the Propagator_DC.m file to create a halo orbit at the Sun-Earth system L2 point with amplitude Az=1,20,000 km. This code uses Richardson's theoretical results as initial conditions for halo orbit propagation using C3BP model and then corrects those initial conditions using the Differential Correction technique to get an accuracy of ~10^-16.

Result
----

HALO ORBIT
----
![Final halo orbit DC=50, tolerances, 2 full circles](https://user-images.githubusercontent.com/61064585/164998765-d52eb5b4-3c57-41ff-86d2-a49a10844fda.png)


Velocity vs Time
----
![Final halo orbit DC=50, tolerances, 2 full circles, velocity profile](https://user-images.githubusercontent.com/61064585/164998770-4658522f-79e0-41df-8c3f-3a78c1fbe793.png)


Chaotic dynamics in the presence of three bodies
----
![deviating considerably from L1, x0=0 788](https://user-images.githubusercontent.com/61064585/164998783-d877e4f6-9095-466a-9044-889b9e9fdbc9.png)
----
Entering into an unstable manifold and becoming heliocentric eventually
