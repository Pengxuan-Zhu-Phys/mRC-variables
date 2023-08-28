# mRC-variables

A set of variables to solve the masses of particles in decay tree with invisible particles. We defined the mRC-variables with the help of the particle physics MC analysis tool: Rivet-3.1.X.

# Documentation

Documentation can be found in [Inspire-HEP](https://inspirehep.net/literature/2182463) and in [arXiv:2211.08132](https://arxiv.org/abs/2211.08132). 

# Demos 
## Installation of Rivet 
The easiest way to install Rivet is via the Docker system, and you can also  install it via the standard script. For example, you can install it in Ubuntu system.  
```bash
sudo apt install texlive-full imagemagick 

cd /scratch/rivet
wget https://gitlab.com/hepcedar/rivetbootstrap/raw/3.1.6/rivet-bootstrap
chmod +x rivet-bootstrap
```  
More information, please see the Rivet homepage in [GitLab](https://gitlab.com/hepcedar/rivet/-/blob/release-3-1-x/doc/tutorials/installation.md). 

## Build analysis 
1. Firstly, you should make a new rivet analysis in your work directory
    ```bash
    rivet-mkanalysis CEPC_SMUON
    ``` 
2. Then you find rivet create three files:
    * `CEPC_SMUON.cc `
    * `CEPC_SMUON.plot`
    * `CEPC_SMUON.info`
3. Download our code and replace the above three files. 
4. Recompile the analysis
    ```bash
    rivet-build RivetCEPC_SMUON.so CEPE_SMUON.cc
    ```

## Analysis the MC events 
The `*.hepmc` event can be analysis in a `.gz` compressed manner. Using command like 

```bash
rivet --pwd -a CEPC_SMUON <YourEventsDir>.hepmc.gz
```

In above command, `—pwd` means Rivet set the environment in current path, `-a CPEC_SMUON` means using the analysis named `CEPC_SMUON`, `<YourEventsDir>.hepmc.gz` is the event file. 
Except the screen output and the output files information defined by yourself, Rivet will output a `Rivet.yoda` file for plotting. Alternatively, you can define the output file name via command like:
```bash
rivet --pwd -o rivet.yoda -a CEPC_SMUON <YourEventsDir>.hepmc.gz
```

## Plot the histograms

```bash
rivet-mkhtml --no-weight Rivet.yoda
```

After this, Rivet output a `rivet-plots` folder which contains a `index.html` and a folder named `CEPC_SMUON`. The histograms figures are stored in the `CEPC_SMUON`.

# 