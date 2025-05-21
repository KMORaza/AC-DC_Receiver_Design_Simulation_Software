using System;

namespace ACDCSimulator
{
    public class PowerFactorCorrector
    {
        private double inputVoltage;
        private double loadResistance;
        private string rectifierType;
        private double frequency;

        public PowerFactorCorrector(double inputVoltage, double loadResistance, string rectifierType, double frequency)
        {
            this.inputVoltage = inputVoltage;
            this.loadResistance = loadResistance;
            this.rectifierType = rectifierType;
            this.frequency = frequency;
        }

        public double CalculateUncorrectedPowerFactor()
        {
            // Simplified power factor calculation based on rectifier type and load
            double phaseAngle = 0;
            if (rectifierType == "Full-Wave")
            {
                // Full-wave rectifier: less phase shift, but non-linear load reduces PF
                phaseAngle = Math.PI / 6; // 30 degrees due to diode conduction
            }
            else if (rectifierType == "Half-Wave")
            {
                // Half-wave: more distortion, worse PF
                phaseAngle = Math.PI / 4; // 45 degrees
            }
            else
            {
                // Default: assume full-wave behavior
                phaseAngle = Math.PI / 6;
            }

            // Adjust phase angle based on load (simplified model)
            double loadEffect = 1 - (100 / (loadResistance + 100)); // Heavy loads reduce PF
            phaseAngle *= (1 + loadEffect * 0.2);

            // Power factor = cos(phase angle)
            double powerFactor = Math.Cos(phaseAngle);
            return Math.Round(powerFactor, 3);
        }

        public double ApplyPowerFactorCorrection()
        {
            // Simulate active PFC (e.g., boost converter)
            double uncorrectedPF = CalculateUncorrectedPowerFactor();
            if (uncorrectedPF >= 0.99)
            {
                // Already near unity, no significant correction needed
                return uncorrectedPF;
            }

            // Active PFC model: boost converter aligns current with voltage
            // Simplified: PFC improves PF to 0.95–0.99 based on input conditions
            double correctionFactor = 0.95 + (0.04 * (inputVoltage / 120)); // Scales with input voltage
            double correctedPF = Math.Min(0.99, uncorrectedPF + (1 - uncorrectedPF) * correctionFactor);

            return Math.Round(correctedPF, 3);
        }

        public double GetCorrectedCurrent(double apparentPower)
        {
            // Calculate corrected current after PFC
            double correctedPF = ApplyPowerFactorCorrection();
            double realPower = apparentPower * correctedPF;
            double correctedCurrent = realPower / inputVoltage;
            return Math.Round(correctedCurrent, 3);
        }
    }
}