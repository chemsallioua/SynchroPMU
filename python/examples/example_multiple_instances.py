from pmu_estimator import PMUEstimator, EstimatorConfig
import math

def create_pmu_estimator(amp, ph, freq):
    sample_rate = 25600
    window_size = 2048
    dt = 1 / sample_rate
    
    # Generating input signal window
    input_signal_window = [amp * math.cos(2 * math.pi * freq * dt * i + ph) for i in range(window_size)]

    # Initializing pmu estimator from config struct
    pmu_config = EstimatorConfig(
        n_cycles = 4,
        f0 = 50,
        frame_rate = 50,
        fs = 25600,
        n_bins = 11,
        P = 3,
        Q = 3,
        interf_trig = 0.0033,
        rocof_thresh = [3.0, 25.0, 0.035],
        rocof_low_pass_coeffs = [0.5913, 0.2043, 0.2043],
        iter_eipdft = True
    )

    # Initializing pmu estimator object
    pmu = PMUEstimator()

    # Configuring pmu estimator from config struct
    if(pmu.configure_from_class(pmu_config) != 0):
        print("Error: Could not configure pmu estimator")
        return None, None

    # setting a mid window fracsec
    mid_window_fracsec = 0.0

    # Estimating frame 
    estimated_frame = pmu.estimate(input_signal_window, mid_window_fracsec)
    if estimated_frame is None:
        print("Error: Could not estimate frame")
        return None, None

    return pmu, estimated_frame

if __name__ == "__main__":

    # Create two PMU estimators with different signal parameters
    pmu1, estimated_frame1 = create_pmu_estimator(2.0, 1.0, 50.0)
    pmu2, estimated_frame2 = create_pmu_estimator(1.5, 0.5, 60.0)

    # Printing estimated frames
    if estimated_frame1 is not None:
        print("Estimated frame from PMU 1:", estimated_frame1)
    
    if estimated_frame2 is not None:
        print("Estimated frame from PMU 2:", estimated_frame2)
    
