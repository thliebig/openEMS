import numpy
from pathlib import Path
import CSXCAD
from openEMS import openEMS
import pandas
import plotly.express as px
from logging import getLogger
from scipy.fft import fft, fftfreq
from scipy.constants import c

logger = getLogger(__name__)

def theoretical_resonance_frequencies(
	# Compute the resonance frequency of a specific mode. See [Wikipedia](https://en.wikipedia.org/wiki/Microwave_cavity#Rectangular_cavity).
	cavity_size_meters:numpy.ndarray,
		# Size of the cavity in meters. E.g. `numpy.array((1,2,3))`.
	mode_number:numpy.ndarray,
		# Mode number. E.g. `numpy.array((1,0,0))` is the fundamental mode in `x` direction.
	E_field_direction:str,
		# Direction of the electric field, either `'x'`, `'y'` or `'z'`.
):
	if sum(mode_number) in {0,1}:
		# This means there is no excitation of the EM field at all, see the expressions of the fields [here](https://www.brainkart.com/article/Rectangular-and-circular-cavity-resonators_12513/).
		return float('NaN')
	match E_field_direction:
		case 'x':
			if mode_number[1]==0 or mode_number[2]==0:
				return float('NaN')
		case 'y':
			if mode_number[0]==0 or mode_number[2]==0:
				return float('NaN')
		case 'z':
			if mode_number[1]==0 or mode_number[0]==0:
				return float('NaN')
		case _:
			raise ValueError(f'Wrong `E_field_direction`. ')
	return c/2*(sum((mode_number/cavity_size_meters)**2))**.5

def configure_rectangular_cavity_simulation(
	# Create an openEMS simulation of a rectangular resonant cavity excited by a Dirac delta to study the resonance modes.
	path_to_folder_where_to_write_output:Path,
		# Folder where to save the results.
	cavity_size:tuple,
		# Size of the cavity, in meters, e.g. `(1,2,3)`.
	f_max:float = 1e9,
		# Highest frequency to be used in the excitation, in Hz.
	frequency_resolution:float = 1e6,
		# Frequency resolution in Hz. Higher resolution require more time steps.
):
	SIMULATION_VOLUME_START = numpy.array(cavity_size)/2 # Position of bottom left cornrer.
	SIMULATION_VOLUME_SIZE = numpy.array(cavity_size) # Simulation volume size.

	# Initialize the machinery:
	CSX = CSXCAD.ContinuousStructure()
	simulation = openEMS()
	simulation.SetCSX(CSX)

	simulation.SetBoundaryCond(['PEC']*6) # "Perfect electric conductor" in the 6 faces of the rectangle, thus creating a resonant cavity.

	simulation.SetDiracExcite(f_max=f_max)
	simulation.SetMaxTime(frequency_resolution**-1)

	# Grid lines:
	CSX.grid.SetDeltaUnit(unit=1) # Set unit to meter.
	for i,(xyz,n_lines) in enumerate(zip(['x','y','z'], [111,111,111])):
		CSX.grid.AddLine(
			xyz,
			list(numpy.linspace(SIMULATION_VOLUME_START[i], SIMULATION_VOLUME_START[i] + SIMULATION_VOLUME_SIZE[i], n_lines)),
		)

	# Excitation port:
	port_size = SIMULATION_VOLUME_SIZE*.1 # Make it small so the excitation is like a point source inside the volume.
	port_in_start = SIMULATION_VOLUME_START + (SIMULATION_VOLUME_SIZE-port_size)/3 # Move the port outside the center of the volume in order to excite a larger number of oscillation modes.
	simulation.AddLumpedPort(
		port_nr = 1,
		R = 1e99, # Infinite resistance
		start = port_in_start,
		stop = (port_in_start + port_size),
		p_dir = 'x',
		excite = 1, # V m⁻¹
	)

	path_to_csx_model_to_save = path_to_folder_where_to_write_output/'CSX_model.xml'
	CSX.Write2XML(path_to_csx_model_to_save)

	logger.info(f'Saved CSX model in: {path_to_csx_model_to_save}')

	return simulation

def run_simulation(simulation, path_to_folder_where_to_write_output:Path):
	try:
		logger.info(f'Starting simulation in {path_to_folder_where_to_write_output}')
		simulation.Run(
			sim_path = path_to_folder_where_to_write_output/'data',
			cleanup = True,
		)
	finally:
		logger.info(f'Simulation finished, results available in {path_to_folder_where_to_write_output}')

def analyze_results(path_to_folder_where_to_write_output:Path, cavity_size:tuple, f_max:float):
	logger.info('Analyzing results...')
	time_domain_data = pandas.read_csv(
		path_to_folder_where_to_write_output/f'data/port_ut_1',
		comment = '%',
		names = ['Time (s)', 'Voltage (V)'],
		sep = '\t',
	)

	MAX_POINTS_TO_PLOT = int(10e3) # To avoid heavy plots, which is a weak point of plotly.
	fig = px.line(
		title = 'Time domain simulation results',
		data_frame = time_domain_data if len(time_domain_data)<MAX_POINTS_TO_PLOT else time_domain_data.sample(n=MAX_POINTS_TO_PLOT).sort_index(),
		x = 'Time (s)',
		y = 'Voltage (V)',
	)
	fig.write_html(path_to_folder_where_to_write_output/'voltage vs time.html')

	frequency_domain_data = pandas.DataFrame(
		{
			'fft': fft(time_domain_data['Voltage (V)'])[0:len(time_domain_data)//2],
			'Frequency (Hz)': fftfreq(len(time_domain_data), numpy.diff(time_domain_data['Time (s)'])[0])[:len(time_domain_data)//2],
		},
	)
	frequency_domain_data = frequency_domain_data.sort_values('Frequency (Hz)')
	frequency_domain_data['abs(fft)'] = numpy.abs(frequency_domain_data['fft'])

	frequency_domain_data = frequency_domain_data.query(f'`Frequency (Hz)` < {f_max}')
	frequency_domain_data = frequency_domain_data.query(f'`Frequency (Hz)` > 0') # Drop the DC component.
	fig = px.line(
		title = 'Frequency domain simulation results<br><sup>Dashed lines indicate theoretical resonance frequencies</sup>',
		data_frame = frequency_domain_data,
		x = 'Frequency (Hz)',
		y = 'abs(fft)',
		log_y = True,
		markers = True,
	)
	fig.update_layout(
		xaxis = dict(
			exponentformat = 'SI'
		)
	)
	# Draw theoretical resonance frequencies:
	for nx in range(5):
		for ny in range(5):
			for nz in range(5):
				resonance_frequency = theoretical_resonance_frequencies(numpy.array(cavity_size),numpy.array((nx,ny,nz)), E_field_direction='x')
				if resonance_frequency > f_max or numpy.isnan(resonance_frequency):
					continue
				fig.add_vline(
					x = resonance_frequency,
					line_dash = 'dot',
				)
	fig.write_html(path_to_folder_where_to_write_output/'frequency_domain.html')

	logger.info(f'Plots available in {path_to_folder_where_to_write_output}')

def main():
	PATH_TO_FOLDER_WHERE_TO_WRITE_SIMULATION_OUTPUT = Path(__file__).parent/'simulation_results'
	CAVITY_SIZE = (1,1,1)
	F_MAX = 1e9

	PATH_TO_FOLDER_WHERE_TO_WRITE_SIMULATION_OUTPUT.mkdir(exist_ok=True, parents=True)

	simulation = configure_rectangular_cavity_simulation(
		path_to_folder_where_to_write_output = PATH_TO_FOLDER_WHERE_TO_WRITE_SIMULATION_OUTPUT,
		cavity_size = CAVITY_SIZE,
		f_max = F_MAX,
	)
	run_simulation(
		simulation = simulation,
		path_to_folder_where_to_write_output = PATH_TO_FOLDER_WHERE_TO_WRITE_SIMULATION_OUTPUT,
	)
	analyze_results(
		path_to_folder_where_to_write_output = PATH_TO_FOLDER_WHERE_TO_WRITE_SIMULATION_OUTPUT,
		cavity_size = CAVITY_SIZE,
		f_max = F_MAX,
	)

if __name__ == '__main__':
	from logging import basicConfig, INFO
	basicConfig(level=INFO, format="%(name)s %(levelname)s: %(message)s")

	main()
