Introduction
************
astroMUGS-ui does exactly what astroMUGS does, i.e. it provides a full post-process pipeline to couple radiative transfer and gas-grain simulations
for models using discretized multiple dust populations. 

Unlike astroMUGS, where the user has to create an instance of their model to create the pipeline, astroMUGS-ui is a node-based interface where the user
can manipulate nodes to build their pipeline intuitively and quickly. The GUI is based on React Flow for the frontend, and fastAPI for the backend. 
The backend directly communicates with the up-to-date astroMUGS package, so that the possibilities are exactly the same as using astroMUGS. 


