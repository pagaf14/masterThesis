# Setting PATH for Python 3.7
# The original version is saved in .zprofile.pysave
#PATH="/Library/Frameworks/Python.framework/Versions/3.7/bin:${PATH}"
#export PATH

# SU2
export PATH=/Users/matteopaganelli/Desktop/masterThesis/SU2/Codes/SU2_bin/bin:$PATH

#SU2 environment variables
export SU2_RUN=/Users/matteopaganelli/Desktop/masterThesis/SU2/Codes/SU2_bin/bin
export SU2_HOME=/Users/matteopaganelli/Desktop/masterThesis/SU2/Codes/SU2
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$PYTHONPATH:$SU2_RUN

# PARAVIEW
export PATH=/Applications/ParaView-5.11.2.app/Contents/MacOS:$PATH
export PATH=/Applications/ParaView-5.11.2.app/Contents/bin:$PATH
export PYTHONPATH="/Applications/ParaView-5.11.2.app/Contents/MacOS:$PYTHONPATH"
export PYTHONPATH="/Applications/ParaView-5.11.2.app/Contents/bin:$PYTHONPATH"

# PY_ENV
eval "$(pyenv init --path)"

# PYCHARM
export TERM=xterm-256color

# MATLAB

export PATH=/Applications/MATLAB_R2024a.app/bin:$PATH

# Setting PATH for Python 3.12
# The original version is saved in .zprofile.pysave
PATH="/Library/Frameworks/Python.framework/Versions/3.12/bin:${PATH}"
export PATH

