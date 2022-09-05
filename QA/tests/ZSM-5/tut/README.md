# Mini tutorial on embedded chemical system

This is a mini tutorial on how to set up calculations for chemical systems
that are embedded in a larger system. The zeolite case that we are interested
in is an example of such a case. Here I show how to do such calculations
but I use an (over simplified) example. 

The example I use here is a chain of Carbon atoms terminated at either with
hydrogen atoms. The Carbon-Carbon bonds are alternating single and triple
bonds. Ultimately we are interested in a reaction that happens at one triple
bond, and we assume that the reaction is constrained by the rest of the 
chemical system (in this example this isn't really true, but for the zeolite
system it is). The example proceeds with three calculations discussed below.

## Input1: optimize the initial full structure

In `input1.nw` I specify a chain of 6 Carbon atoms, terminated at either end
with a Hydrogen atom. This structure is optimized. Comparing this to our
zeolite challenge `input1.nw` is the step that gives us the actual full
zeolite structure. For the zeolite case we don't have to do this step 
because because we have been given the structure already.

## Input2: truncate and terminate the model structure

In `input2.nw` we remove the end CH groups, and replace the new Carbon ends
with Hydrogens (named h_c because these hydrogens replace Carbons). This 
way we obtain a model that contains just the central triple-bonded C-C 
group. Because the newly introduced Hydrogen atoms are not at equilibrium
positions we need to optimize those. However, we want to keep the central 
Carbon atoms in their current positions as they are supposed to model
Carbon atoms in a long chain. By setting
```
   set geometry:actlist 1 4
```
we tell NWChem to optimize only the positions of the terminal Hydrogen
atoms.

## Input3: the system reacts subject to a constraining environment

In `input3.nw` we have a singlet Oxygen atom react with the Carbon-Carbon
triple bond. In this calculation the terminal Hydrogen atoms model the
rest of the original Carbon chain that constrain the movement of the
reacting Carbon atoms. Setting
```
   set geometry:actlist 2 3 5
```
allows the central Carbon atoms to move as well as the Oxygen atom, but 
the Hydrogen atoms are kept fixed.

## Our zeolite challenge

For our zeolite system we essentially need to generate inputs 2 and 3. 
So first we need to carve out a cluster of the zeolite around the active
site. The active site is of course the AlOH piece, and we need to include
some of the silicate around that. Then we need to terminate this cluster
with hydrogens, and optimize the positions of the Hydrogens. 

In the next step we will keep the terminal Hydrogens fixed to model the
effect of the larger zeolite structure. We will introduce the propanol
and model the reaction.

Please let me know if you have any questions.

