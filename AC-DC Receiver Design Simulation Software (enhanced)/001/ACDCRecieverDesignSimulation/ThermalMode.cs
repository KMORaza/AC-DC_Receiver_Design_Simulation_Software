using System;

namespace ACDCSimulator
{
    public class ThermalModel
    {
        private double junctionTemperature; // in °C
        private double ambientTemperature; // in °C
        private double thermalResistance; // Junction-to-ambient in °C/W
        private double thermalCapacity; // in J/°C
        private double lastUpdateTime;

        public double JunctionTemperature => junctionTemperature;
        public double PowerDissipation { get; private set; }

        public ThermalModel(double ambientTemp = 25, double thermalResistance = 10, double thermalCapacity = 0.1)
        {
            ambientTemperature = ambientTemp;
            this.thermalResistance = thermalResistance;
            this.thermalCapacity = thermalCapacity;
            junctionTemperature = ambientTemp;
            lastUpdateTime = 0;
        }

        public void Update(double powerDissipation, double currentTime)
        {
            PowerDissipation = powerDissipation;

            // Simple thermal model: RC network approximation
            double deltaTime = currentTime - lastUpdateTime;
            if (deltaTime <= 0) return;

            // Calculate temperature change using simplified thermal model
            double steadyStateTemp = ambientTemperature + (powerDissipation * thermalResistance);
            double timeConstant = thermalResistance * thermalCapacity;

            // Exponential approach to steady state temperature
            junctionTemperature = steadyStateTemp - (steadyStateTemp - junctionTemperature) *
                                Math.Exp(-deltaTime / timeConstant);

            lastUpdateTime = currentTime;
        }

        public double CalculateDegradationFactor()
        {
            // Calculate performance degradation based on temperature
            const double maxTemp = 150; // Maximum operating temperature in °C
            const double refTemp = 25; // Reference temperature in °C

            if (junctionTemperature <= refTemp)
                return 1.0;

            double tempRatio = (junctionTemperature - refTemp) / (maxTemp - refTemp);
            return Math.Max(0.5, 1.0 - 0.5 * tempRatio); // Linear degradation from 100% to 50%
        }
    }
}