This directory contains MATLAB scripts for modeling gastric motility and gut-brain axis interactions. Below is a description of each file:

--------------------------------------------------
1. compartment_gut_brain_axis.m
   - Models the vagal region of the gut-brain axis influence on gastric motility and gastric emptying.
   - Incorporates neural control mechanisms (vagal) affecting gastric activity.

2. compartment_radius.m
   - Simulates compartmental radius under deformation

3. diff_sol.m
   - Solves differential equations governing stomach dynamics.
   - Uses numerical solvers (e.g., `ode23s`, `ode15s`) for system integration.

--------------------------------------------------
# USAGE INSTRUCTIONS:
- Run `diff_sol.m` to simulate the stomach model.
- If you want gastric volume to be constant, comment out Volume_total variable in 'compartment_gut_brain_axis.m'
