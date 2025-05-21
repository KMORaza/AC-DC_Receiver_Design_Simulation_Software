using System;

namespace ACDCSimulator
{
    public class ThermalModel
    {
        private double junctionTemperature; // in °C
        private double caseTemperature; // in °C
        private double heatSinkTemperature; // in °C
        private double ambientTemperature; // in °C
        private double thermalResistanceJunctionToCase; // °C/W
        private double thermalResistanceCaseToHeatSink; // °C/W
        private double thermalResistanceHeatSinkToAmbient; // °C/W
        private double thermalCapacityJunction; // J/°C
        private double thermalCapacityHeatSink; // J/°C
        private double lastUpdateTime;
        private bool thermalRunawayDetected;

        public double JunctionTemperature => junctionTemperature;
        public double CaseTemperature => caseTemperature;
        public double HeatSinkTemperature => heatSinkTemperature;
        public double PowerDissipation { get; private set; }
        public bool ThermalRunawayDetected => thermalRunawayDetected;

        public ThermalModel(
            double ambientTemp = 25,
            double thermalResistanceJunctionToCase = 1.0,
            double thermalResistanceCaseToHeatSink = 0.5,
            double thermalResistanceHeatSinkToAmbient = 10.0,
            double thermalCapacityJunction = 0.01,
            double thermalCapacityHeatSink = 10.0)
        {
            ambientTemperature = ambientTemp;
            this.thermalResistanceJunctionToCase = thermalResistanceJunctionToCase;
            this.thermalResistanceCaseToHeatSink = thermalResistanceCaseToHeatSink;
            this.thermalResistanceHeatSinkToAmbient = thermalResistanceHeatSinkToAmbient;
            this.thermalCapacityJunction = thermalCapacityJunction;
            this.thermalCapacityHeatSink = thermalCapacityHeatSink;
            junctionTemperature = ambientTemp;
            caseTemperature = ambientTemp;
            heatSinkTemperature = ambientTemp;
            lastUpdateTime = 0;
            thermalRunawayDetected = false;
        }

        public void Update(double powerDissipation, double currentTime)
        {
            PowerDissipation = powerDissipation;
            double deltaTime = currentTime - lastUpdateTime;
            if (deltaTime <= 0) return;

            // Update heat sink temperature
            double heatSinkSteadyState = ambientTemperature + powerDissipation * thermalResistanceHeatSinkToAmbient;
            double heatSinkTimeConstant = thermalResistanceHeatSinkToAmbient * thermalCapacityHeatSink;
            heatSinkTemperature = heatSinkSteadyState - (heatSinkSteadyState - heatSinkTemperature) *
                                Math.Exp(-deltaTime / heatSinkTimeConstant);

            // Update case temperature
            double caseSteadyState = heatSinkTemperature + powerDissipation * thermalResistanceCaseToHeatSink;
            caseTemperature = caseSteadyState - (caseSteadyState - caseTemperature) *
                              Math.Exp(-deltaTime / (thermalResistanceCaseToHeatSink * thermalCapacityHeatSink));

            // Update junction temperature
            double junctionSteadyState = caseTemperature + powerDissipation * thermalResistanceJunctionToCase;
            double junctionTimeConstant = thermalResistanceJunctionToCase * thermalCapacityJunction;
            junctionTemperature = junctionSteadyState - (junctionSteadyState - junctionTemperature) *
                                  Math.Exp(-deltaTime / junctionTimeConstant);

            // Check for thermal runaway (threshold: 150°C)
            thermalRunawayDetected = junctionTemperature > 150;

            lastUpdateTime = currentTime;
        }

        public double CalculateDegradationFactor()
        {
            const double maxTemp = 150; // Maximum operating temperature in °C
            const double refTemp = 25; // Reference temperature in °C

            if (junctionTemperature <= refTemp)
                return 1.0;

            double tempRatio = (junctionTemperature - refTemp) / (maxTemp - refTemp);
            return Math.Max(0.5, 1.0 - 0.5 * tempRatio);
        }
    }
}