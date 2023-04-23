# OpticalTweezerProject

The script calculates the significant force fluctuations caused by breaking or shifting of actin bundles, and produces a histogram of the force differences found in the data. 

<!-- GETTING STARTED -->
## Getting Started

The script requires Python 3, and the following libraries:

- pandas
- numpy
- scipy
- matplotlib

### Installation

_Below is an instruction of how you can install and run the code._

1. Clone the repo
   ```sh
   git clone https://github.com/JulianMeigen/OpticalTweezerProject.git
   cd OpticalTweezerProject
   ```
2. Install Python packages
   ```sh
   pip install -r requirements.txt
   ```
3. Run the python script
   ```sh
   python OT_Project.py -f 20220919/processed_curves-20220919 # requires the path for your data
   ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage
Here are a few examples how you can use the script:
1. Calculates the significant force fluctuations caused by breaking or shifting of actin bundles for a singular file and saves a plot.
   ```sh
   python OT_Project.py -f 20220919/processed_curves-20220919/force-save-2022.09.19-14.09.55.589 # requires the path for your datafile
   ```
2. Calculates the force fluctuations for all datafiles (.txt) in the given folder and produces a histogram containing all force differences.
   ```sh
   python OT_Project.py -f 20220919/processed_curves-20220919 # requires the path for your data
   ```
3. Use different coulumns as x and y axis
   ```sh
   python OT_Project.py -f [FILENAME] -x "distance" -y "xSignal1" # default paramaters
   python OT_Project.py -f 20220919/processed_curves-20220919 -x "time" -y "xSignal2"
   ```
4. Use your own Calibration data of the Optical Tweezer
   ```sh
   python OT_Project.py -f 20220919/processed_curves-20220919 -c TRUE
   ```
   
<p align="right">(<a href="#readme-top">back to top</a>)</p>

