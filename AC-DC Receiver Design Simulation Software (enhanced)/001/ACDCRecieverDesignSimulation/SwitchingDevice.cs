using System;

namespace ACDCSimulator
{
    public enum SwitchingDeviceType { MOSFET, IGBT, GaN, SiC }

    public abstract class SwitchingDevice
    {
        public double Voltage { get; protected set; }
        public double Current { get; protected set; }
        public double OnResistance { get; protected set; }
        public double SwitchingFrequency { get; protected set; }
        public double GateCharge { get; protected set; }
        public ThermalModel ThermalModel { get; protected set; }
        public double Efficiency { get; protected set; }

        public abstract void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt);
        public abstract double CalculateConductionLoss(double loadCurrent);
        public abstract double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency);

        protected void InitializeThermalModel(double ambientTemp = 25)
        {
            ThermalModel = new ThermalModel(ambientTemp);
        }

        protected void UpdateThermalModel(double totalLoss, double currentTime)
        {
            ThermalModel.Update(totalLoss, currentTime);
            Efficiency = 1.0 - (totalLoss / (Voltage * Current + 1e-9)); // Avoid division by zero
        }
    }

    public class MOSFET : SwitchingDevice
    {
        public MOSFET(double ambientTemp = 25)
        {
            OnResistance = 0.05; // ohms
            GateCharge = 50e-9; // Coulombs
            SwitchingFrequency = 100e3; // Hz
            InitializeThermalModel(ambientTemp);
        }

        public override void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt)
        {
            double time = DateTime.Now.Ticks / 1e7; // Current time in seconds

            // Calculate losses
            double conductionLoss = CalculateConductionLoss(loadCurrent);
            double switchingLoss = CalculateSwitchingLoss(loadCurrent, inputVoltage, SwitchingFrequency);
            double totalLoss = conductionLoss + switchingLoss;

            // Update thermal model
            UpdateThermalModel(totalLoss, time);

            // Apply temperature-dependent degradation
            double degradationFactor = ThermalModel.CalculateDegradationFactor();
            double effectiveOnResistance = OnResistance / degradationFactor;

            Voltage = inputVoltage * dutyCycle;
            Current = loadCurrent;
            Voltage -= Current * effectiveOnResistance; // Account for on-state voltage drop
        }

        public override double CalculateConductionLoss(double loadCurrent)
        {
            return loadCurrent * loadCurrent * OnResistance;
        }

        public override double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency)
        {
            return 0.5 * inputVoltage * loadCurrent * GateCharge * frequency;
        }
    }

    public class IGBT : SwitchingDevice
    {
        public IGBT()
        {
            OnResistance = 0.1; // ohms
            GateCharge = 100e-9; // Coulombs
            SwitchingFrequency = 50e3; // Hz
        }

        public override void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt)
        {
            Voltage = inputVoltage * dutyCycle;
            Current = loadCurrent;
            Voltage -= 1.5 + Current * OnResistance; // V_CE(sat) + resistive drop
        }

        public override double CalculateConductionLoss(double loadCurrent)
        {
            return 1.5 * loadCurrent + loadCurrent * loadCurrent * OnResistance;
        }

        public override double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency)
        {
            return 0.75 * inputVoltage * loadCurrent * GateCharge * frequency;
        }
    }

    public class GaN : SwitchingDevice
    {
        public GaN()
        {
            OnResistance = 0.02; // ohms
            GateCharge = 20e-9; // Coulombs
            SwitchingFrequency = 500e3; // Hz
        }

        public override void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt)
        {
            Voltage = inputVoltage * dutyCycle;
            Current = loadCurrent;
            Voltage -= Current * OnResistance;
        }

        public override double CalculateConductionLoss(double loadCurrent)
        {
            return loadCurrent * loadCurrent * OnResistance;
        }

        public override double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency)
        {
            return 0.3 * inputVoltage * loadCurrent * GateCharge * frequency;
        }
    }

    public class SiC : SwitchingDevice
    {
        public SiC()
        {
            OnResistance = 0.03; // ohms
            GateCharge = 30e-9; // Coulombs
            SwitchingFrequency = 200e3; // Hz
        }

        public override void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt)
        {
            Voltage = inputVoltage * dutyCycle;
            Current = loadCurrent;
            Voltage -= Current * OnResistance;
        }

        public override double CalculateConductionLoss(double loadCurrent)
        {
            return loadCurrent * loadCurrent * OnResistance;
        }

        public override double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency)
        {
            return 0.4 * inputVoltage * loadCurrent * GateCharge * frequency;
        }
    }
}