from pmu_estimator import PMUEstimator, EstimatorConfig
import math

if __name__ == "__main__":

    # Input Signal parameters
    AMP = 2.0
    PH = 1.0
    FREQ = 50.0
    sample_rate = 25600
    window_size = 2048
    dt = 1 / sample_rate
    
    # Generating input signal window
    input_signal_window = [AMP * math.cos(2 * math.pi * FREQ * dt * i + PH) for i in range(window_size)]

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
    pmu.configure_from_class(pmu_config)

    # setting a mid window fracsec
    mid_window_fracsec = 0.0

    # Estimating frame
    estimated_frame = pmu.estimate(input_signal_window, mid_window_fracsec)

    # Printing estimated frame
    print(estimated_frame)
