# ARLas
Auditory Research Lab audio software (ARLas)

ABOUT ARLas

ARLas is a MATLAB-based software package written by Shawn Goodman, director of the Auditory Research Laboratory in the Department of Communication Sciences and Disorders at the University of Iowa. It was primarily designed to record otoacoustic emissions, but may be suitable for use with any application requiring synchronous averaging (e.g. for auditory evoked potential recording). The software is very flexible and allows the user to customize experiment protocols. The software does not provide any data analysis. It is assumed that the user is familiar with MATLAB and can write his or her own analysis code and apply it post hoc. 

DISCLAIMER

THE SOFTWARE IS PROVIDED TO THE RESEARCH COMMUNITY "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

COPYRIGHT

Copyright © 2011-2025 Shawn S. Goodman. If you have not been given this software by the author, please email him at 
shawn-goodman@uiowa.edu to let him know that you are using the software.  Please also be so kind as to acknowledge the use of this software in any publications that include data acquired using this software.

The ARLas software uses a piece of third-party software, called playrec.  Playrec is a MATLAB and Octave utility (MEX file) that provides access to soundcards using PortAudio, a free, open-source audio I/O library. Copyright © 2006-2008 Robert Humphrey.  Permission is hereby granted, free of charge, to any person obtaining a copy of this software [playrec] and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  Playrec may be accessed online at http://www.playrec.co.uk/index.php.


SYSTEM  REQUIREMENTS

Operating System
ARLas was written for and tested on both Windows operating systems and Macintosh operating systems. It will run on both 32- and 64-bit operating systems. In theory, it should also be useable on Linux, but this has not been explicitly tested.  

MATLAB
ARLas is written using object oriented programming.  It was originally developed on and tested on MATLAB version 24.1 (R2024a).  Some slightly older versions may work, but much older versions may not, due to changes in the way MATLAB handles classes and objects.  The MATLAB signal processing toolbox must be installed in order to run ARLas. Both 32- and 64-bit versions of MATLAB will support ARLas.

Sound Card
ARLas should be able to access any soundcard, including on-board cards, PCI-based cards, and external sound cards (connected using USB or Firewire).  Soundcards can be accessed through the playrec utility via different host API including ASIO, WMME and DirectSound (under Windows) and Core Audio (Macintosh).  

While ARLas will run using any card, the ability to record otoacoustic emissions requires a high-quality 24-bit card and up-to-date drivers.  ARLas works best using ASIO drivers, if they are available.  One of the limitations of ARLas arises from the third-party software, playrec, which it uses.  Playrec does not have the ability to specify the number of bits.  In rare cases, it has been observed that certain drivers can cause a card with variable bit rate capability, including 24 bits, to run in 16-bit mode.  It is critical that this is avoided, because 16 bits does not provide a low enough noise floor for recording otoacoustic emissions.  As of this writing, ASIO drivers have not caused this problem, and therefore they are preferred.  
