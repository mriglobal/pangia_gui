# Non-Dockerized Instructions:

# INITIAL INSTALL:
# In terminal cd to this directory.
# Run "conda env create -f gui.yml".
# Run "flask db init".
# Run "flask db migrate –m".
# Run "flask db upgrade".
# Run "conda install -c bioconda fastp".

# STARTING THE SERVICE:
# In three separate terminals cd to this directory.
# In each terminal run "conda activate gui".
# In 1st terminal run "redis-server".
# In 2nd terminal run "rq worker pangia-tasks".
# In 3rd terminal run "export FLASK_APP=pangia_gui.py".
#    			- then "flask run".

# ACCESSING/USING THE SERVICE:
# With terminals running, navigate to "localhost:5000".

# DEBUGING
# Go to "http://localhost:5000/auth/reset_site" if you need to reset the database.

# Dockerized Instructions:

# INITIAL INSTALL:
# Follow all Non-Dockerized INITIAL INSTALL instructions.
# Locate and open 'docker-compose.yml' within this directory.
# Locate the PATH specified under: services --> app --> volumes. 
	Change this line such that it corresponds with your PanGIA directory. Save and close the file.
# Run "docker-compose build".

# STARTING THE SERVICE:
# In terminal cd to this directory.
# Run "docker-compose up".
# With compose container running, navigate to "localhost:8000".