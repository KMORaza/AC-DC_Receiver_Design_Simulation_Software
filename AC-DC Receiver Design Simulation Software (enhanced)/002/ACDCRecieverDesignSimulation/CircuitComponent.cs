using MathNet.Numerics.LinearAlgebra;
using System;
using System.Numerics;

namespace ACDCSimulator
{
    public enum RectifierType { HalfWave, FullWave, Bridge }
    public enum FilterType { None, Capacitive, Inductive, Active }
    public enum RegulatorType { None, Linear, Switching }

    public class Diode
    {
        public double ForwardVoltage { get; set; } = 0.7;
        public double ReverseBreakdown { get; set; } = 100;
        public double Current { get; private set; }
        public double Voltage { get; private set; }

        public void Update(double voltage)
        {
            Voltage = voltage;
            Current = voltage > ForwardVoltage ? Math.Exp((voltage - ForwardVoltage) / 0.025) : 0;
            if (voltage < -ReverseBreakdown) Current = -1e-6;
        }
    }

    public class Capacitor
    {
        public double Capacitance { get; set; }
        public double Voltage { get; private set; }
        public double Current { get; private set; }

        public void Update(double current, double dt)
        {
            Current = current;
            Voltage += current * dt / Capacitance;
        }
    }

    public class Inductor
    {
        public double Inductance { get; set; }
        public double Current { get; private set; }
        public double Voltage { get; private set; }

        public void Update(double voltage, double dt)
        {
            Voltage = voltage;
            Current += voltage * dt / Inductance;
        }
    }

    public class Transformer
    {
        public double TurnsRatio { get; set; }
        public double PrimaryVoltage { get; set; }
        public double SecondaryVoltage { get; private set; }

        public void Update()
        {
            SecondaryVoltage = PrimaryVoltage * TurnsRatio;
        }
    }

    public class ZenerDiode
    {
        public double ZenerVoltage { get; set; }
        public double ForwardVoltage { get; set; } = 0.7;
        public double Voltage { get; private set; }
        public double Current { get; private set; }

        public void Update(double voltage)
        {
            Voltage = voltage;
            if (voltage > ForwardVoltage)
                Current = Math.Exp((voltage - ForwardVoltage) / 0.025);
            else if (voltage < -ZenerVoltage)
                Current = -Math.Exp((-voltage - ZenerVoltage) / 0.025);
            else
                Current = 0;
        }
    }

    public class OpAmp
    {
        public double Gain { get; set; } = 10000;
        public double OutputVoltage { get; private set; }

        public void Update(double inputPlus, double inputMinus)
        {
            OutputVoltage = Gain * (inputPlus - inputMinus);
            OutputVoltage = Math.Max(Math.Min(OutputVoltage, 15), -15);
        }
    }

    public class PIDController
    {
        public double Kp { get; set; } // Proportional gain
        public double Ki { get; set; } // Integral gain
        public double Kd { get; set; } // Derivative gain
        private double integral = 0;
        private double previousError = 0;

        public double Update(double error, double dt)
        {
            integral += error * dt;
            double derivative = (error - previousError) / dt;
            previousError = error;
            return Kp * error + Ki * integral + Kd * derivative;
        }
    }

    public class PLL
    {
        public double ReferenceFrequency { get; set; }
        public double Kp { get; set; } = 0.1; // Phase detector gain
        public double Ki { get; set; } = 0.01; // Integrator gain
        private double phaseErrorIntegral = 0;
        private double vcoPhase = 0;

        public double Update(double inputPhase, double dt)
        {
            double phaseError = inputPhase - vcoPhase;
            phaseErrorIntegral += phaseError * dt;
            double vcoFrequency = ReferenceFrequency + Kp * phaseError + Ki * phaseErrorIntegral;
            vcoPhase += 2 * Math.PI * vcoFrequency * dt;
            return vcoPhase;
        }
    }

    public class PWM
    {
        public double CarrierFrequency { get; set; } = 1000; // Hz
        public double DutyCycle { get; private set; }

        public void Update(double controlSignal, double maxVoltage)
        {
            // Simplified triangular wave comparison
            double carrier = Math.Sin(2 * Math.PI * CarrierFrequency * DateTime.Now.Ticks / 1e7);
            DutyCycle = controlSignal / maxVoltage;
            DutyCycle = Math.Max(0, Math.Min(1, DutyCycle));
        }
    }

    public class DigitalControlLoop
    {
        public double SamplingRate { get; set; } = 1000; // Hz
        public double[] Coefficients { get; set; } // For digital filter (e.g., FIR)
        private double[] previousInputs;
        private double[] previousOutputs;

        public DigitalControlLoop(int order)
        {
            previousInputs = new double[order];
            previousOutputs = new double[order];
            Coefficients = new double[order];
        }

        public double Update(double input)
        {
            // Shift previous inputs/outputs
            for (int i = previousInputs.Length - 1; i > 0; i--)
            {
                previousInputs[i] = previousInputs[i - 1];
                previousOutputs[i] = previousOutputs[i - 1];
            }
            previousInputs[0] = input;

            // Compute output (simplified FIR filter)
            double output = 0;
            for (int i = 0; i < Coefficients.Length; i++)
            {
                output += Coefficients[i] * previousInputs[i];
            }
            previousOutputs[0] = output;
            return output;
        }
    }

    public class CircuitParameters
    {
        public RectifierType Rectifier { get; set; } = RectifierType.Bridge;
        public FilterType Filter { get; set; } = FilterType.None;
        public RegulatorType Regulator { get; set; } = RegulatorType.None;
        public SwitchingDeviceType SwitchingDevice { get; set; } = SwitchingDeviceType.MOSFET;
        public double TransformerTurnsRatio { get; set; } = 1.0;
        public double FilterCapacitance { get; set; } = 100e-6;
        public double FilterInductance { get; set; } = 10e-3;
        public double FilterCutoffFrequency { get; set; } = 100;
        public double RegulatorVoltage { get; set; } = 5.0;
        public double PID_Kp { get; set; } = 1.0;
        public double PID_Ki { get; set; } = 0.1;
        public double PID_Kd { get; set; } = 0.01;
        public double PLL_ReferenceFrequency { get; set; } = 60;
        public double PWM_CarrierFrequency { get; set; } = 1000;
        public double[] DigitalCoefficients { get; set; } = new double[] { 1.0 };
    }

    public class SimulationResult
    {
        public double SteadyStateDcVoltage { get; set; }
        public double RippleVoltage { get; set; }
        public double PowerDissipation { get; set; }
        public double DominantHarmonicFrequency { get; set; }
        public double DominantHarmonicAmplitude { get; set; }
        public double TotalHarmonicDistortion { get; set; } // Added for THD
        public double[] TransientTimePoints { get; set; }
        public double[] TransientVoltagePoints { get; set; }
        public double[] BodeFrequencyPoints { get; set; }
        public double[] BodeMagnitudePoints { get; set; }
        public double[] BodePhasePoints { get; set; }
        public double[] NyquistRealPoints { get; set; }
        public double[] NyquistImagPoints { get; set; }
        public Complex[] RootLocusPoints { get; set; }
        public double DeviceTemperature { get; set; }
        public double DeviceEfficiency { get; set; }
    }
}