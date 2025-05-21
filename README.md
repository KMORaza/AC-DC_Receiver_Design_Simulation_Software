## AC/DC Receiver Design Simulation Software

_**(May 20 2025)**_ This software provides an environment for simulation and analysis of AC/DC receiver design. It has a wide range of functionalities and it works fine but I honestly say that it isn't as effective as I wanted it to be.

_**(May 21 2025)**_ I have integrated a few improvements in this software and I have written the [updated codebase](https://github.com/KMORaza/AC-DC_Receiver_Design_Simulation_Software/tree/main/AC-DC%20Receiver%20Design%20Simulation%20Software%20(enhanced)) so it is better than before.

---

The software simulates an AC-to-DC power conversion system, allowing users to configure circuit parameters (rectifier, filter, regulator, switching device, etc.), run transient simulations, and analyze results through metrics (output voltage, ripple, power dissipation) and plots (waveform, Bode, Nyquist, root locus). It models components like diodes, capacitors, inductors, transformers, and control systems (PID, PLL, PWM) using numerical methods and provides thermal, impedance, noise, and power factor analyses.

1. **Purpose of the Software**: The AC/DC Receiver Simulator is a Windows Forms application designed to simulate and analyze the performance of an AC-to-DC power conversion circuit, allowing users to configure components like rectifiers, filters, regulators, and control systems, and to evaluate metrics such as output voltage, ripple, power dissipation, harmonic distortion, impedance, noise, and stability through transient and frequency-domain analyses.

2. **User Interface Initialization (Form1.cs)**: The application initializes a Windows Forms interface with two tabs—Simulation and Analysis—containing input fields for circuit parameters (e.g., input voltage, load resistance, rectifier type), a "Run Simulation" button, result labels (e.g., output voltage, ripple), and charts for visualizing transient waveforms, Bode plots, Nyquist plots, and root locus plots; dropdowns for rectifier, filter, regulator, and switching device types are populated with predefined enums, and charts are cleared to prevent initialization errors.

3. **Input Validation (Form1.cs)**: When the "Run Simulation" button is clicked, the software validates user inputs, ensuring that numerical values (e.g., input voltage, load resistance, turns ratio, PID gains) are positive and properly formatted, and that digital coefficients are comma-separated valid numbers; invalid inputs trigger error messages, halting the simulation.

4. **Circuit Configuration (CircuitParameters Class)**: The software collects user inputs into a `CircuitParameters` object, which stores settings such as rectifier type (HalfWave, FullWave, Bridge), filter type (None, Capacitive, Inductive, Active), regulator type (None, Linear, Switching), switching device type (MOSFET, IGBT, GaN, SiC), transformer turns ratio, filter capacitance (default 100 µF), filter inductance (default 10 mH), filter cutoff frequency (default 100 Hz), regulator voltage (default 5 V), PID gains (Kp=1, Ki=0.1, Kd=0.01), PLL reference frequency (default 60 Hz), PWM carrier frequency (default 1000 Hz), and digital filter coefficients (default [1.0]).

5. **Simulation Core (CircuitSimulator.cs)**: The `CircuitSimulator` class is instantiated with the input voltage, load resistance, and `CircuitParameters` object, initializing simulation parameters including a fixed input frequency of 60 Hz, a time step (dt) of 10 µs, and 10,000 simulation steps, corresponding to a total simulation time of 0.1 seconds.

6. **Component Initialization**: The simulator initializes circuit components based on user selections: a `Transformer` with the specified turns ratio, an array of `Diode` objects (1 for HalfWave, 2 for FullWave, 4 for Bridge), a `Capacitor` for capacitive filters, an `Inductor` for inductive filters, an `OpAmp` for active filters, a `ZenerDiode` for linear regulators, a `PIDController` with user-defined gains, a `PLL` with the reference frequency, a `PWM` with the carrier frequency, a `DigitalControlLoop` with user-defined coefficients, and a `SwitchingDevice` (MOSFET, IGBT, GaN, or SiC) for switching regulators; initial device temperature is set to 25°C, and efficiency to 100%.

7. **Transient Simulation Loop**: The simulation runs for 10,000 steps, with each step calculating the system state at time t = i * dt (i from 0 to 9999); the loop updates component states, computes voltages and currents, and accumulates power dissipation, storing time and output voltage points for transient analysis.

8. **Transformer Operation**: The `Transformer` sets its primary voltage as V_primary = inputVoltage * sin(2π * 60 * t) and computes the secondary voltage as V_secondary = V_primary * TurnsRatio, modeling an ideal transformer without losses.

9. **Rectifier Operation**: The rectifier processes the transformer’s secondary voltage: for HalfWave, a single diode conducts when V_secondary > 0.7 V, yielding V_rectified = V_secondary - 0.7; for FullWave, two diodes handle positive and negative cycles, yielding V_rectified = max(V_secondary - 0.7, -V_secondary - 0.7); for Bridge, four diodes produce V_rectified = |V_secondary| - 1.4, accounting for two diode drops; diode current is modeled as I = exp((V - 0.7)/0.025) for forward bias or -1 µA for reverse breakdown.

10. **PLL Frequency Tracking**: The `PLL` tracks the input signal’s phase (inputPhase = 2π * 60 * t) by computing phase error (inputPhase - vcoPhase), updating the phase error integral, and adjusting the VCO frequency as f_vco = ReferenceFrequency + Kp * phaseError + Ki * phaseErrorIntegral (Kp=0.1, Ki=0.01); the VCO phase is updated as vcoPhase += 2π * f_vco * dt, and the rectified voltage is scaled by cos(pllPhase - inputPhase) to account for phase alignment.

11. **Filter Operation (Capacitive)**: For a capacitive filter, the capacitor current is computed as I = (V_rectified - V_output) / loadResistance, and the capacitor voltage (output voltage) is updated as V_output += I * dt / Capacitance, modeling the charging/discharging dynamics.

12. **Filter Operation (Inductive)**: For an inductive filter, the inductor voltage is V_rectified - V_output, and the inductor current is updated as I += (V_rectified - V_output) * dt / Inductance; the output voltage is V_output = V_rectified - V_inductor, where V_inductor is derived from the inductor’s state.

13. **Filter Operation (Active)**: For an active filter, an `OpAmp` processes the rectified voltage as the non-inverting input, with the output voltage as the inverting input; the op-amp output is V_opamp = Gain * (V_rectified - V_output) (Gain=10,000), clipped between -15 V and 15 V; the output voltage is updated as V_output += (V_opamp - V_output) * (1 - exp(-2π * FilterCutoffFrequency * dt)), approximating a first-order low-pass filter.

14. **PID Control**: The `PIDController` computes the error as error = RegulatorVoltage - V_output, updates the integral as integral += error * dt, and the derivative as derivative = (error - previousError) / dt; the PID output is Kp * error + Ki * integral + Kd * derivative, and the output voltage is adjusted as V_output += PID_output * dt.

15. **PWM Modulation**: The `PWM` computes a duty cycle as DutyCycle = V_output / RegulatorVoltage, clamped between 0 and 1, using a simplified triangular wave comparison with a carrier signal sin(2π * CarrierFrequency * t); the output voltage is scaled by the duty cycle, V_output *= DutyCycle.

16. **Digital Control Loop**: The `DigitalControlLoop` implements a simplified FIR filter, shifting previous inputs and outputs, adding the current input, and computing the output as a weighted sum: output = Σ(Coefficients[i] * previousInputs[i]), updating the output voltage as V_output = output.

17. **Switching Regulator**: For a switching regulator, a `SwitchingDevice` (MOSFET, IGBT, GaN, or SiC) is updated with the PWM duty cycle, input voltage (rectified voltage), load current (V_output / loadResistance), and dt; the device computes conduction and switching losses, updates its thermal model, and adjusts the output voltage as V_output = switchDevice.Voltage, accounting for on-state voltage drops.

18. **Linear Regulator**: For a linear regulator, a `ZenerDiode` limits the output voltage to the Zener voltage (RegulatorVoltage) if V_output exceeds it, modeling ideal regulation; the Zener current is computed as I = exp((V - 0.7)/0.025) for forward bias or -exp((-V - ZenerVoltage)/0.025) for reverse bias.

19. **Switching Device Models (SwitchingDevice.cs)**: Each switching device (MOSFET, IGBT, GaN, SiC) has specific parameters: MOSFET (R_on=0.05 Ω, Q_g=50 nC, f_sw=100 kHz), IGBT (R_on=0.1 Ω, Q_g=100 nC, f_sw=50 kHz, V_CE(sat)=1.5 V), GaN (R_on=0.02 Ω, Q_g=20 nC, f_sw=500 kHz), SiC (R_on=0.03 Ω, Q_g=30 nC, f_sw=200 kHz); each device computes conduction and switching losses and updates its thermal model.

20. **Conduction Loss Calculation**: For MOSFET, GaN, and SiC, conduction loss is I_load² * R_on; for IGBT, it is 1.5 * I_load + I_load² * R_on, accounting for the collector-emitter saturation voltage.

21. **Switching Loss Calculation**: Switching loss is computed as 0.5 * V_in * I_load * Q_g * f_sw for MOSFET, 0.75 * V_in * I_load * Q_g * f_sw for IGBT, 0.3 * V_in * I_load * Q_g * f_sw for GaN, and 0.4 * V_in * I_load * Q_g * f_sw for SiC, reflecting differences in switching characteristics.

22. **Thermal Model (ThermalModel.cs)**: The `ThermalModel` tracks junction, case, and heat sink temperatures using a thermal network: heat sink temperature is updated as T_hs = T_ambient + P * R_th_hs-a - (T_ambient + P * R_th_hs-a - T_hs) * exp(-dt / (R_th_hs-a * C_th_hs)); case temperature as T_case = T_hs + P * R_th_c-hs - (T_hs + P * R_th_c-hs - T_case) * exp(-dt / (R_th_c-hs * C_th_hs)); junction temperature as T_j = T_case + P * R_th_j-c - (T_case + P * R_th_j-c - T_j) * exp(-dt / (R_th_j-c * C_th_j)); thermal runaway is detected if T_j > 150°C.

23. **Degradation Factor**: The thermal model computes a degradation factor as 1 - 0.5 * (T_j - 25)/(150 - 25) for T_j > 25°C, minimum 0.5, increasing effective on-resistance as R_on / degradationFactor to model temperature-dependent performance.

24. **Efficiency Calculation**: Device efficiency is computed as 1 - (totalLoss / (V * I + 1e-9)), where totalLoss is conductionLoss + switchingLoss, expressed as a percentage in the results.

25. **Power Dissipation**: Total power dissipation is accumulated as powerSum += (V_output² / loadResistance + conductionLoss + switchingLoss) * dt during the simulation, and the average is computed as powerSum / (steps * dt), stored in the `SimulationResult`.

26. **Steady-State and Ripple Analysis**: The steady-state DC voltage is the average of the output voltage over the last 50% of simulation steps (steps/2 to steps-1), and the ripple voltage is the difference between the maximum and minimum voltages in the same period, reflecting the filter’s effectiveness.

27. **Frequency-Domain Analysis (FFT)**: The output voltage array is transformed using a Fast Fourier Transform (FFT) with MathNet.Numerics’ Fourier.Forward, converting voltages to complex numbers; the dominant harmonic is identified by finding the maximum amplitude (excluding DC) in the first half of the spectrum, with frequency f = index / (steps * dt) and amplitude A = magnitude / steps.

28. **Total Harmonic Distortion (THD)**: THD is computed by identifying the fundamental frequency (60 Hz, index = 60 * steps * dt), summing the squared amplitudes of harmonics (2nd to 10th, at indices n * fundamentalIndex) as harmonicSumSquared, and calculating THD = sqrt(harmonicSumSquared) / fundamentalAmplitude * 100, expressed as a percentage.

29. **Frequency Response (Bode Plot)**: The frequency response is computed over 100 logarithmically spaced frequencies from 0.01 Hz to 10 kHz; for each frequency f, the complex frequency s = j * 2π * f is used to calculate the transfer function G(s) = (Kp + Ki/s + Kd*s) * filter, where filter = 1 / (s * C) for capacitive filters or 1 otherwise; magnitude is 20 * log10(|G|), phase is arg(G) * 180/π, stored in `SimulationResult`.

30. **Nyquist Plot**: The Nyquist plot is generated using the real and imaginary parts of G(s) for the same frequency points, plotting G(s).Real vs. G(s).Imaginary to assess system stability.

31. **Root Locus**: The root locus is approximated as -1 / G(s) for each frequency point, providing a simplified representation of pole movement, stored as complex numbers in `SimulationResult`.

32. **Impedance Matching (ImpedanceMatcher.cs)**: The `ImpedanceMatcher` calculates input impedance as loadResistance / (TurnsRatio²) * rectifierFactor (0.5 for HalfWave, 0.8 for FullWave, 1.0 for Bridge) and output impedance as loadResistance * filterFactor * regulatorFactor, where filterFactor accounts for capacitive (1 / (1 + R/X_C)) or inductive (1 + X_L/R) effects at 60 Hz, and regulatorFactor is 0.9 for linear or 0.7 for switching regulators; reflection coefficient is |(Z_in - Z_source)/(Z_in + Z_source)|, and mismatch loss is -10 * log10(1 - Γ²).

33. **Noise Analysis (NoiseAndDistortionAnalyzer.cs)**: The signal-to-noise ratio (SNR) is computed as 10 * log10(signalPower / totalNoisePower), where signalPower = V_output² / loadResistance, and totalNoisePower = k * T * B * rectifierNoiseFactor (k=1.38e-23 J/K, T=300 K, B=filterCutoff, rectifierNoiseFactor=1.2 for FullWave, 1.5 otherwise).

34. **THD in Noise Analysis**: THD is calculated as sqrt((A_2² + A_3²)) / A_fundamental * 100, where A_fundamental = V_output, A_2 = A_fundamental * harmonicFactor (0.05 for FullWave, 0.1 otherwise), and A_3 = A_2 * 0.5, providing an alternative THD estimate.

35. **EMI Calculation**: Electromagnetic interference (EMI) is estimated as baseEMI + 20 * log10(pwmCarrierFreq / 1000) - 10 * log10(1 + L * C * 10⁻⁹ * (2π * filterCutoff)²), where baseEMI = 100 dBµV, reflecting PWM switching effects and filter attenuation.

36. **Power Factor Correction (PowerFactorCorrector.cs)**: The uncorrected power factor is cos(phaseAngle), with phaseAngle = π/6 (30°) for FullWave, π/4 (45°) for HalfWave, adjusted by loadEffect = 1 - (100 / (loadResistance + 100)) * 1.2; active PFC improves the power factor to min(0.99, uncorrectedPF + (1 - uncorrectedPF) * (0.95 + 0.04 * (inputVoltage / 120))).

37. **Thermal Runaway Check**: After simulation, the switching device’s thermal model is checked for thermal runaway (T_j > 150°C); if detected, the status label displays “Simulation Complete (Thermal Warning)” in orange, otherwise “Simulation Complete” in lime green.

38. **Result Display**: The simulation results are displayed in the UI: output voltage (V), ripple voltage (V), power dissipation (W), dominant harmonic (Hz, V), THD (%), power factor, SNR (dB), EMI (dBµV), input impedance (Ω), output impedance (Ω), reflection coefficient, mismatch loss (dB), device temperature (°C), device efficiency (%), and thermal status (Safe or Warning).

39. **Transient Waveform Plot**: The transient waveform is plotted on `chartWaveform` as a line chart, with x-axis as time (t * 1000 ms) and y-axis as output voltage (V), using black color for visibility.

40. **Bode Plot**: The Bode plot is plotted on `chartBode` with two series: magnitude (dB, black line) and phase (degrees, dark gray line) vs. frequency (Hz), using a logarithmic x-axis from 0.1 Hz to 100 kHz; the chart is hidden if no valid frequency points are available.

41. **Nyquist Plot**: The Nyquist plot is plotted on `chartNyquist` as a point chart, with x-axis as real part and y-axis as imaginary part of G(s), using black points to indicate stability characteristics.

42. **Root Locus Plot**: The root locus is plotted on `chartRootLocus` as a point chart, with x-axis as real part and y-axis as imaginary part of -1/G(s), using black points to approximate pole locations.

43. **Status Updates**: The status label (`lblStatus`) updates during simulation: “Simulating…” (white) during execution, “Simulation Complete” (lime green) or “Simulation Complete (Thermal Warning)” (orange) on success, or “Simulation Failed” (red) on exceptions, ensuring user feedback.

44. **Error Handling**: Exceptions during simulation (e.g., numerical errors, invalid configurations) are caught, displaying an error message with the exception details, setting the status to “Simulation Failed,” and re-enabling the simulate button.

45. **Asynchronous Simulation**: The simulation runs asynchronously using `Task.Run` to prevent UI freezing, allowing the interface to remain responsive during computation-intensive tasks.

---

![](https://github.com/KMORaza/AC-DC_Receiver_Design_Simulation_Software/blob/main/AC-DC%20Receiver%20Design%20Simulation%20Software%20(enhanced)/002/ACDCRecieverDesignSimulation/screenshot.png?raw=true)
