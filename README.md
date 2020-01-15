# ParseEVPSC

### Version 3.3. Jan 2020

Fixed bug while reading multi-step processes

### Version 3.2. Jan 2020

The previous version crashed while parsing some calculation results. In fact HP-EVPSC or EVPSC have slightly different output formats. Some versions output how much of the deformation is elastic vs. plastic. Others do not.

The version 3.2 of parseEVPSC should work with both types of output, as well as 3 and 4 indices for miller indices.
