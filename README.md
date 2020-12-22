# mTowerEventSelectionNew
Updated version of mTowerEventSelection. Includes implementation of side-by-side selection using criteria-based cuts (developed by A. Bochove and R. Bosley) and an adapted kT algorithm (developed by H. Yokoyama.)
This code requires a working installation of ROOT, GCC 7 or higher, and a working install of fastjet. It is recommended that this code be run using lxplus, however a guide to adapting for outside lxplus exists below.

#Quickly run Event Selection from lxplus
1. Firstly, you will need to amend lines 14 - 18 of Analyse_mTower_Run1250_combined.sh to point to your local repository of mTowerEventSelectionNew. The two fastjet directories in lines 12 and 13 should be accessible on lxplus, and will not require amending (though you can find a different fastjet installation if you want to.)
2. If you wish to run a run other than Run 1250, you will need to amend line 19 of Analyse_mTower_Run1250_combined.sh to read 'Analyse_mTower(<YourRunNumberHere>) (you may wish to copy this to a more appropriately named shell script.)
3. Simply run 'Analyse_mTower_Run1250_combined.sh':

```bash
./Analyse_mTower_Run1250_combined.sh
```

4. Firstly, a recent build of ILCSoft from CLICdp on cvmfs will initialise. This is not required, but is useful for ensuring you are using the correct versions of ROOT and GCC. To remove this step, remove line 4 of Analyse_mTower_Run1250_combined.sh.
5. You will see various class definitions compile (these are classes we are importing from the 'classes' directory). Finally, Analyse_mTower_combine.cxx will compile, and run. Do not be alarmed if you see several warnings about potentially uninitialised variables here - this is normal.
6. The code will likely take several hours to run. You may wish to run the code via the batch system, in which case follow the directions in the instructions below.
7. You should find two files have been produced - results_Run_1250_C2_10p_C4_120p_C6_50pixelsC246_3layerclusteringfixed.root, and Run_1250_eventsLeft_normalwindow_C246_layer2fixedclusters.root.
8. Run_1250_eventsLeft_normalwindow_C246_layer2fixedclusters.root is the main output file, and contains three trees:
  
  Frames - contains a single branch 'nHits,' which describes the distribution of total hits per event in events which were selected by Robbie's selection.
  
  tree - contains many branches describing the output from Hiroki's selection, including the X and Y of the shower, and the status of the selection.
  
  prepreparedselections - contains three branches:
    
    event - describes the event number
    
    HirokiStatus - an int_t describing whether or not the associated event was selected by Hiroki's algorithm (0 if not selected, 1 if selected)
    
    RobbieStatus - an int_t describing whether or not the associated event was selected by Robbie's algorithm (0 if not selected, 1 if selected)
9. prepreparedselections is the tree that should be used to quickly find which events were selected by the two algorithms.


# Running Event Selection using the lxplus Condor batch system
1. A script already exists in order to run the code via the lxplus batch system. To run it, simply use:

```bash
condor_submit analyse_Run1250_combined.sub
```

2. To check the status of the job, use:

```bash
condor_q
```

3. You may wish to alter line 1 of the steering file analyse_Run1250_combined.sub for a different shell script if you are using a run other than Run 1250, and copy the file to give it a more appropriate name.


# Adapting to running outside lxplus
1. Firstly, ensure that you have a working installation of ROOT, and GCC v7 or above.
2. If you do not have an installation of cvmfs in your shell, you will need to delete line 4 of Analyse_mTower_Run1250_combined.sh, to avoid trying to initialise ILCSoft.
3. Ensure you have a working installation of Fastjet. If you need to install Fastjet, follow the guide at http://fastjet.fr/quickstart.html.
4. Amend lines 12 and 13 of Analyse_mTower_Run1250_combined.sh, to point to your installation of Fastjet. Note: your fastjet-<VersionNumber>/ directory is not applicable, you will need the fastjet-install/ directory.
5. Ensure you have a local copy of H. Yokoyama's utility files. They are stored on EOS under /eos/project/m/mtower/public/hiroki/00_util/. It is recommended to copy them from here using scp or rsynch, e.g.:

```bash
scp -rp <user>@lxplus.cern.ch:/eos/project/m/mtower/public/hiroki/00_util <PathToCopyDirectory>
```

6. You will also need a copy of 'cellinfo.h' from /eos/project/m/mtower/public/hiroki/11_id_pileup/

7. You will then need to amend lines 32 - 39 of Analyse_mTower_combine.cxx to reflect your local copies of the utility files copied from Hiroki's repository on EOS.
8. Of course, you will not be able to run on CERN's batch system while not in the lxplus shell, so you will have to adapt to your local batch system, or run the code interactively.

