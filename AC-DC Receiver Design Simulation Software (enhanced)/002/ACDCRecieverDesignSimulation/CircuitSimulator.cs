using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Linq;
using System.Numerics;

namespace ACDCSimulator
{
    public class CircuitSimulator
    {
        private readonly double inputVoltage;
        private readonly double loadResistance;
        private readonly CircuitParameters parameters;
        private readonly double frequency = 60;
        private readonly double dt = 1e-5;
        private readonly int steps = 10000;

        public CircuitSimulator(double inputVoltage, double loadResistance, CircuitParameters parameters)
        {
            this.inputVoltage = inputVoltage;
            this.loadResistance = loadResistance;
            this.parameters = parameters;
        }

        public SimulationResult RunSimulation()
        {
            var result = new SimulationResult();
            var timePoints = new double[steps];
            var voltagePoints = new double[steps];

            // Initialize components
            var transformer = new Transformer { TurnsRatio = parameters.TransformerTurnsRatio };
            var diodes = parameters.Rectifier == RectifierType.Bridge ? new Diode[4] :
                         parameters.Rectifier == RectifierType.FullWave ? new Diode[2] : new Diode[1];
            for (int i = 0; i < diodes.Length; i++) diodes[i] = new Diode();
            var capacitor = parameters.Filter == FilterType.Capacitive ? new Capacitor { Capacitance = parameters.FilterCapacitance } : null;
            var inductor = parameters.Filter == FilterType.Inductive ? new Inductor { Inductance = parameters.FilterInductance } : null;
            var opAmp = parameters.Filter == FilterType.Active ? new OpAmp() : null;
            var zener = parameters.Regulator == RegulatorType.Linear ? new ZenerDiode { ZenerVoltage = parameters.RegulatorVoltage } : null;
            var pid = new PIDController { Kp = parameters.PID_Kp, Ki = parameters.PID_Ki, Kd = parameters.PID_Kd };
            var pll = new PLL { ReferenceFrequency = parameters.PLL_ReferenceFrequency };
            var pwm = new PWM { CarrierFrequency = parameters.PWM_CarrierFrequency };
            var digitalLoop = new DigitalControlLoop(parameters.DigitalCoefficients.Length) { Coefficients = parameters.DigitalCoefficients };
            result.DeviceTemperature = 25; // Initial ambient temperature
            result.DeviceEfficiency = 100; // Initial 100% efficiency

            // Initialize switching device for switching regulator
            SwitchingDevice switchDevice = null;
            if (parameters.Regulator == RegulatorType.Switching)
            {
                switch (parameters.SwitchingDevice)
                {
                    case SwitchingDeviceType.MOSFET:
                        switchDevice = new MOSFET();
                        break;
                    case SwitchingDeviceType.IGBT:
                        switchDevice = new IGBT();
                        break;
                    case SwitchingDeviceType.GaN:
                        switchDevice = new GaN();
                        break;
                    case SwitchingDeviceType.SiC:
                        switchDevice = new SiC();
                        break;
                }
            }

            double outputVoltage = 0;
            double powerSum = 0;
            double[] voltagesForRipple = new double[steps];

            // Transient simulation
            for (int i = 0; i < steps; i++)
            {
                double t = i * dt;
                timePoints[i] = t;

                // Input AC signal
                transformer.PrimaryVoltage = inputVoltage * Math.Sin(2 * Math.PI * frequency * t);
                transformer.Update();
                double rectifiedVoltage = 0;

                // Rectifier
                switch (parameters.Rectifier)
                {
                    case RectifierType.HalfWave:
                        diodes[0].Update(transformer.SecondaryVoltage);
                        rectifiedVoltage = diodes[0].Current > 0 ? transformer.SecondaryVoltage - diodes[0].ForwardVoltage : 0;
                        break;
                    case RectifierType.FullWave:
                        diodes[0].Update(transformer.SecondaryVoltage);
                        diodes[1].Update(-transformer.SecondaryVoltage);
                        rectifiedVoltage = Math.Max(diodes[0].Current > 0 ? transformer.SecondaryVoltage - diodes[0].ForwardVoltage : 0,
                                                  diodes[1].Current > 0 ? -transformer.SecondaryVoltage - diodes[1].ForwardVoltage : 0);
                        break;
                    case RectifierType.Bridge:
                        diodes[0].Update(transformer.SecondaryVoltage);
                        diodes[1].Update(-transformer.SecondaryVoltage);
                        diodes[2].Update(transformer.SecondaryVoltage);
                        diodes[3].Update(-transformer.SecondaryVoltage);
                        rectifiedVoltage = Math.Abs(transformer.SecondaryVoltage) - 2 * diodes[0].ForwardVoltage;
                        break;
                }

                // PLL for frequency tracking
                double inputPhase = 2 * Math.PI * frequency * t;
                double pllPhase = pll.Update(inputPhase, dt);
                rectifiedVoltage *= Math.Cos(pllPhase - inputPhase); // Adjust voltage based on phase alignment

                // Filter
                if (capacitor != null)
                {
                    double current = (rectifiedVoltage - outputVoltage) / loadResistance;
                    capacitor.Update(current, dt);
                    outputVoltage = capacitor.Voltage;
                }
                else if (inductor != null)
                {
                    inductor.Update(rectifiedVoltage - outputVoltage, dt);
                    outputVoltage = rectifiedVoltage - inductor.Voltage;
                }
                else if (opAmp != null)
                {
                    double inputVoltage = rectifiedVoltage;
                    opAmp.Update(inputVoltage, outputVoltage);
                    outputVoltage += (opAmp.OutputVoltage - outputVoltage) * (1 - Math.Exp(-2 * Math.PI * parameters.FilterCutoffFrequency * dt));
                }
                else
                {
                    outputVoltage = rectifiedVoltage;
                }

                // PID control
                double error = parameters.RegulatorVoltage - outputVoltage;
                double pidOutput = pid.Update(error, dt);
                outputVoltage += pidOutput * dt;

                // PWM modulation
                pwm.Update(outputVoltage, parameters.RegulatorVoltage);
                outputVoltage *= pwm.DutyCycle;

                // Digital control loop
                outputVoltage = digitalLoop.Update(outputVoltage);

                // Switching device for switching regulator
                if (switchDevice != null)
                {
                    double loadCurrent = outputVoltage / loadResistance;
                    switchDevice.Update(pwm.DutyCycle, rectifiedVoltage, loadCurrent, dt);
                    outputVoltage = switchDevice.Voltage;
                    double conductionLoss = switchDevice.CalculateConductionLoss(loadCurrent);
                    double switchingLoss = switchDevice.CalculateSwitchingLoss(loadCurrent, rectifiedVoltage, parameters.PWM_CarrierFrequency);
                    powerSum += (conductionLoss + switchingLoss) * dt;
                }
                // Regulator (linear or fallback for switching without device)
                else if (zener != null)
                {
                    zener.Update(outputVoltage);
                    outputVoltage = Math.Min(outputVoltage, zener.ZenerVoltage);
                }
                else if (parameters.Regulator == RegulatorType.Switching)
                {
                    double dutyCycle = parameters.RegulatorVoltage / outputVoltage;
                    dutyCycle = Math.Max(0, Math.Min(1, dutyCycle));
                    outputVoltage *= dutyCycle;
                }

                voltagePoints[i] = outputVoltage;
                voltagesForRipple[i] = outputVoltage;
                powerSum += outputVoltage * outputVoltage / loadResistance * dt;
            }

            // Store thermal results after simulation loop
            if (switchDevice != null)
            {
                result.DeviceTemperature = switchDevice.ThermalModel.JunctionTemperature;
                result.DeviceEfficiency = switchDevice.Efficiency * 100; // Convert to percentage
            }

            // Analyze results
            result.TransientTimePoints = timePoints;
            result.TransientVoltagePoints = voltagePoints;
            result.SteadyStateDcVoltage = voltagesForRipple.Skip(steps / 2).Average();
            result.RippleVoltage = voltagesForRipple.Skip(steps / 2).Max() - voltagesForRipple.Skip(steps / 2).Min();
            result.PowerDissipation = powerSum / (steps * dt);

            // Frequency-domain analysis (FFT)
            var complexVoltages = voltagesForRipple.Select(v => new Complex32((float)v, 0)).ToArray();
            Fourier.Forward(complexVoltages, FourierOptions.Matlab);
            double maxAmplitude = 0;
            int maxIndex = 0;
            for (int i = 1; i < complexVoltages.Length / 2; i++)
            {
                double amplitude = complexVoltages[i].Magnitude;
                if (amplitude > maxAmplitude)
                {
                    maxAmplitude = amplitude;
                    maxIndex = i;
                }
            }
            result.DominantHarmonicFrequency = maxIndex / (steps * dt);
            result.DominantHarmonicAmplitude = maxAmplitude / steps;

            // Total Harmonic Distortion (THD) calculation
            double fundamentalFreq = frequency; // 60 Hz
            int fundamentalIndex = (int)(fundamentalFreq * steps * dt);
            double fundamentalAmplitude = fundamentalIndex < complexVoltages.Length ? complexVoltages[fundamentalIndex].Magnitude / steps : 0;
            double harmonicSumSquared = 0;
            const int maxHarmonicOrder = 10; // Consider up to 10th harmonic
            for (int n = 2; n <= maxHarmonicOrder; n++)
            {
                int harmonicIndex = fundamentalIndex * n;
                if (harmonicIndex < complexVoltages.Length / 2)
                {
                    double harmonicAmplitude = complexVoltages[harmonicIndex].Magnitude / steps;
                    harmonicSumSquared += harmonicAmplitude * harmonicAmplitude;
                }
            }
            result.TotalHarmonicDistortion = fundamentalAmplitude > 0 ? Math.Sqrt(harmonicSumSquared) / fundamentalAmplitude * 100 : 0;

            // Stability and frequency response analysis
            ComputeFrequencyResponse(result);

            return result;
        }

        private void ComputeFrequencyResponse(SimulationResult result)
        {
            // Simplified transfer function for the system (e.g., PID + filter)
            int numPoints = 100;
            double[] frequencies = new double[numPoints];
            double[] magnitudes = new double[numPoints];
            double[] phases = new double[numPoints];
            double[] nyquistReal = new double[numPoints];
            double[] nyquistImag = new double[numPoints];
            Complex[] rootLocus = new Complex[numPoints];

            double fMin = 1e-2, fMax = 1e4;
            for (int i = 0; i < numPoints; i++)
            {
                double f = fMin * Math.Pow(fMax / fMin, (double)i / (numPoints - 1));
                frequencies[i] = f;
                Complex s = new Complex(0, 2 * Math.PI * f);

                // Simplified transfer function: G(s) = PID + Filter
                Complex pid = parameters.PID_Kp + parameters.PID_Ki / s + parameters.PID_Kd * s;
                Complex filter = parameters.Filter == FilterType.Capacitive ? 1 / (s * parameters.FilterCapacitance) : 1;
                Complex G = pid * filter;

                magnitudes[i] = 20 * Math.Log10(G.Magnitude);
                phases[i] = G.Phase * 180 / Math.PI;
                nyquistReal[i] = G.Real;
                nyquistImag[i] = G.Imaginary;

                // Root locus (simplified: poles of 1 + k*G(s) = 0)
                rootLocus[i] = -1 / G; // Placeholder for root locus
            }

            result.BodeFrequencyPoints = frequencies;
            result.BodeMagnitudePoints = magnitudes;
            result.BodePhasePoints = phases;
            result.NyquistRealPoints = nyquistReal;
            result.NyquistImagPoints = nyquistImag;
            result.RootLocusPoints = rootLocus;
        }
    }
}