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

        public abstract void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt);
        public abstract double CalculateConductionLoss(double loadCurrent);
        public abstract double CalculateSwitchingLoss(double loadCurrent, double inputVoltage, double frequency);
    }

    public class MOSFET : SwitchingDevice
    {
        public MOSFET()
        {
            OnResistance = 0.05; // ohms
            GateCharge = 50e-9; // Coulombs
            SwitchingFrequency = 100e3; // Hz
        }

        public override void Update(double dutyCycle, double inputVoltage, double loadCurrent, double dt)
        {
            Voltage = inputVoltage * dutyCycle;
            Current = loadCurrent;
            Voltage -= Current * OnResistance; // Account for on-state voltage drop
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