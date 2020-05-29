# Neuron Receptive Field Finding

Neurons in most of biological visual system responds to stimuli in local areas in the visual field, this area is termed receptive field by neuroscientist. Similar phenomenon has been described for artificial neural networks like CNN. [See Scholarpedia](http://www.scholarpedia.org/article/Receptive_field)

To find receptive field, neuroscientists usually display images on different part of the visual field, and record the respond of the neuron. This process is repeated several times to reduce neuronal noise. This repo is dedicated to infer the position of receptive field from an array of such noisy recordings. 

## Basic Plan
To conduct statistical inference on the receptive field, we would like to have a forward model of Stimulus-response relationship of visual neurons and the noise distribution of neuronal firing. 

The most basic model of visual neuron is just a linear filter $K$ with a nonlinearity of firing threshold $\phi$. So first we assume the filter is a gaussian kernel and disregard the threshold. 

As for noise, neuronal firing process could be modelled as a Poisson process with varying underlying rate. In large rate approximation, the noise will be Gaussian like. With smaller rate we can use Poisson model to evaluate the likelihood of each observation. 

As long as we have the likelihood, we can optimize the likelihood function with various methods. One particular method implemented in this repo is MCMC, i.e. to sample from the posterior distribution of the right parameters. Then we see the distribution of the parameters. This process could also be thought of as a kind of stochastic optimization. 
