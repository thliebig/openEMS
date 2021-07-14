Installing `openEMS` with `Octave` in `Windows 10`:

1. Be aware that `Octave` does not allow spaces in path. For example, users shoul d not install `Octave` nor should they not install `openEMS` into `Program Files` folder. Since this path name contains a ' ' white space, the user will encounter problem trying to access `openEMS` from within `Octave`. 

2. When both `openEMS` and `Octave` are present on the `Windows` computer, the user will need to an extra step to 'link' `openEMS` up with `Octave`. This step uses `Octave` startup file to tell `Octave` where to locate `openEMS`. The user will need to find two instances of `octaverc` file. 

  - The files are located like so [see "Octave Support" section](http://openems.de/index.php/OpenEMS#Installation):

    - `Octave install directory`\share\octave\site\m\startup\octaverc
    - `Octave install directory`\share\octave\4.0.0\m\startup\octaverc

  - On a Windows machine, user should be aware that `mingw` or `mingw64` is used, in which case the paths to find these files will be:

    - `Octave install directory`\ __mingw64__ \share\octave\site\m\startup\octaverc
    - `Octave install directory`\ __mingw64__ \share\octave\5.2.0\m\startup\octaverc

  3. The user will need to manually edit the `octaverc` files adding this line. For convenience, the line can be appended to the end of the file. 
  
  ```
  addpath(" `openEMS executable location` /matlab/");
  ```
  
  4. When the above steps are completed, follow the steps outlined [in this URL](http://openems.de/index.php/Tutorial:_First_Steps).
  
