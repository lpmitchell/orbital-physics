using Avionic.Math;

namespace Avionic.Simulation.Kepler
{
    public class CelestialBody
    {
        public readonly double Mass;
        public readonly OrbitalElements OrbitalElements;

        public CelestialBody(double mass, OrbitalElements orbit)
        {
            Mass = mass;
            OrbitalElements = orbit;
        }
    }
}