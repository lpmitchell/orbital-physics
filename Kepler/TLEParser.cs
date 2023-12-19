using System;
using System.Globalization;
using Avionic.Math;
using UnityEngine;

namespace Avionic.Simulation.Kepler
{
    public class TLEParser
    {
        public static OrbitalData ToOrbitalData(CelestialBody parent, string tleLine1, string tleLine2)
        {
            var epochYear = int.Parse(tleLine1.Substring(18, 2));
            var epochDay = double.Parse(tleLine1.Substring(20, 12), CultureInfo.InvariantCulture);
            var inclination = double.Parse(tleLine2.Substring(8, 8), CultureInfo.InvariantCulture);
            var longitudeOfAscendingNode = double.Parse(tleLine2.Substring(17, 8), CultureInfo.InvariantCulture);
            var eccentricity = double.Parse("0." + tleLine2.Substring(26, 7), CultureInfo.InvariantCulture);
            var argumentOfPeriapsis = double.Parse(tleLine2.Substring(34, 8), CultureInfo.InvariantCulture);
            var meanAnomaly = double.Parse(tleLine2.Substring(43, 8), CultureInfo.InvariantCulture);
            var meanMotion = double.Parse(tleLine2.Substring(52, 11), CultureInfo.InvariantCulture);
            
            var meanMotionPerSecond = (24 * 60 * 60) / meanMotion;
            
            var epoch = ConvertYearDayToDateTime(epochYear, epochDay);
            
            var gravitationalParameter = MathConstants.G * parent.Mass;
            var semiMajorAxis = System.Math.Cbrt(System.Math.Pow((meanMotionPerSecond / (2 * System.Math.PI)), 2) * gravitationalParameter);

            var correctedMeanAnomaly = meanAnomaly * System.Math.Sqrt(gravitationalParameter / System.Math.Pow(semiMajorAxis, 3));
            correctedMeanAnomaly %= 2 * System.Math.PI;
            
            var orbit = new OrbitalData(
                eccentricity, 
                semiMajorAxis, 
                correctedMeanAnomaly, 
                inclination, 
                argumentOfPeriapsis,
                longitudeOfAscendingNode, 
                parent.Mass).CloneToUnityCoordinateSpace();

            orbit.Epoch = epoch;
            
            return orbit;
        }
        
        private static DateTime ConvertYearDayToDateTime(int year, double dayOfYear)
        {
            var baseYear = (year < 57) ? 2000 : 1900; // TLE format uses two-digit years
            var dateTime = new DateTime(baseYear + year, 1, 1).AddDays(dayOfYear - 1);
            return dateTime;
        }

    }
}