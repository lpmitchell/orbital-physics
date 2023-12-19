using Avionic.Math;
using static System.Math;

namespace Avionic.Simulation.Kepler
{
    public struct OrbitalElements
    {
        // Gravitational constant
        private const double G = 6.67430e-11;
        
        // Unit vector
        private static readonly Vector3d K = new(0, 0, 1);

        // Used for comparisons
        private const double ParabolaAccuracy = 1e-9;
        private const double Accuracy = 1e-9;

        public readonly CelestialBody Parent;
        public readonly Vector3d Position;
        public readonly Vector3d Velocity;

        private readonly double _positionMag;
        private readonly double _velocityMag;
        private readonly double _trueAnomaly;
        private readonly Vector3d _ascendingNodeVector;
        private readonly double _longitudeOfAscendingNode;
        private readonly double _argumentOfPeriapsis;
        private readonly double _gravitationalParameter;
        private readonly double _specificMechanicalEnergy;
        private readonly Vector3d _angularMomentumVector;
        private readonly double _angularMomentum;
        private readonly Vector3d _eccentricityVector;
        private readonly double _eccentricAnomaly;
        private readonly double _eccentricity;
        private readonly double _alpha;
        private readonly double _semiMajorAxis;
        private readonly double _orbitalPeriod;
        private readonly double _meanMotion;
        private readonly double _semilatusRectum;
        private readonly double _inclination;
        private readonly double _meanAnomaly;

        private readonly bool _isParabolic;
        private readonly bool _isHyperbolic;
        private readonly bool _isElliptical;

        public bool isParabolic => _isParabolic;
        public bool isElliptical => _isElliptical;
        public bool isHyperbolic => _isHyperbolic;
        public bool isCircular => Abs(_eccentricity) < 0.00001;
        
        public OrbitalElements(CelestialBody parent, Vector3d position, Vector3d velocity)
        {
            Parent = parent;
            Position = position;
            Velocity = velocity;
            
            // Populate struct here!
            _positionMag = Position.magnitude;
            _velocityMag = Velocity.magnitude;
            
            // Calculations ---
            _gravitationalParameter = G * Parent.Mass;
            _angularMomentumVector = Vector3d.Cross(Position, Velocity);
            _angularMomentum = _angularMomentumVector.magnitude;
            
            _eccentricityVector = CalculateEccentricityVector(_gravitationalParameter, _velocityMag, _positionMag, Position, Velocity);
            _eccentricity = _eccentricityVector.magnitude;
            if (_eccentricity < Accuracy) _eccentricity = 0; 
            
            
            _ascendingNodeVector = CalculateAscendingNodeVector(_angularMomentumVector);
            _longitudeOfAscendingNode = CalculateLongitudeOfAscendingNode(_ascendingNodeVector);
            _trueAnomaly = CalculateTrueAnomaly(_eccentricityVector, _eccentricity, Position, Velocity, _positionMag);
            _eccentricAnomaly = CalculateEccentricAnomaly(_eccentricity, _trueAnomaly);
            _argumentOfPeriapsis = CalculateArgumentOfPeriapsis(_ascendingNodeVector, _eccentricityVector, _eccentricity);
            _specificMechanicalEnergy = (Pow(_velocityMag, 2) / 2) - (_gravitationalParameter / _positionMag);
            
            var alpha = -_specificMechanicalEnergy * 2d / _gravitationalParameter;
            if (Abs(alpha) < ParabolaAccuracy)
            {
                _alpha = 0;
            }
            else
            {
                _alpha = alpha;
            }
            
            _isParabolic = Abs(_alpha) < Accuracy;
            _isHyperbolic = _alpha < 0 && !_isParabolic;
            _isElliptical = _alpha > 0 && !_isParabolic;

            _meanAnomaly = CalculateMeanAnomaly(_isElliptical, _isHyperbolic, _eccentricAnomaly, _eccentricity);
            _semiMajorAxis = CalculateSemiMajorAxis(_isParabolic, _alpha);
            _orbitalPeriod = CalculateOrbitalPeriod(_isElliptical, _semiMajorAxis, _gravitationalParameter);
            _meanMotion = 2 * PI / _orbitalPeriod;
            
            _semilatusRectum = _semiMajorAxis * (1 - Pow(_eccentricity, 2));
            _inclination = Acos(_angularMomentumVector.z / _angularMomentumVector.magnitude);
        }

        public double OrbitalPeriod => _orbitalPeriod;

        public double GetApoapsis()
        {
            return _semiMajorAxis * (1 + _eccentricity);
        }
        
        public double GetPeriapsis()
        {
            return _semiMajorAxis * (1 - _eccentricity);
        }
        
        public double DeltaTimeToPeriapsis()
        {
            var dT = _orbitalPeriod / (2 * PI) * -_meanAnomaly; // M_p is 0
            return dT >= 0 ? dT : dT + _orbitalPeriod;
        }

        public double DeltaTimeToApoapsis()
        {
            var dT = DeltaTimeToPeriapsis() + _orbitalPeriod / 2;
            if (dT > _orbitalPeriod) dT -= _orbitalPeriod;
            return dT;
        }


        private static double CalculateTrueAnomaly(Vector3d eccentricityVector, double eccentricity, Vector3d position, Vector3d velocity, double positionMag)
        {
            if (Abs(eccentricity) < Accuracy)
            {
                return 0;
            }
            
            var v = Acos(Vector3d.Dot(eccentricityVector, position) / (eccentricity * positionMag));
            if (Vector3d.Dot(position, velocity) < 0)
            {
                v = 2 * PI - v;
            }

            return v;
        }

        private static double CalculateMeanAnomaly(bool isElliptical, bool isHyperbolic, double eccentricAnomaly, double eccentricity)
        {
            if (!isElliptical && !isHyperbolic) return 0;
            
            if (isElliptical)
            {
                // For elliptical orbits
                return eccentricAnomaly - eccentricity * Sin(eccentricAnomaly);
            }

            
            // For hyperbolic orbits
            return eccentricity * Sinh(eccentricAnomaly) - eccentricAnomaly;
        }


        private static Vector3d CalculateAscendingNodeVector(Vector3d angularMomentumVector)
        {
            var n = Vector3d.Cross(K, angularMomentumVector);
            return n.magnitude < Accuracy ? new Vector3d(1, 0, 0) : n; 
        }

        private static double CalculateLongitudeOfAscendingNode(Vector3d ascendingNodeVector)
        {
            if (ascendingNodeVector.magnitude <= Accuracy)
            {
                return 0;
            }
            
            var longitude = Acos(ascendingNodeVector.x / ascendingNodeVector.magnitude);
            if (ascendingNodeVector.y < 0)
            {
                longitude = 2 * PI - longitude;
            }
            return longitude;
        }


        private static double CalculateArgumentOfPeriapsis(Vector3d ascendingNodeVector, Vector3d eccentricityVector, double eccentricity)
        {
            var omega = Acos(Vector3d.Dot(ascendingNodeVector, eccentricityVector) / (ascendingNodeVector.magnitude * eccentricity));
            if (eccentricityVector.z < 0)
            {
                omega = 2 * PI - omega;
            }

            return omega;
        }
        
        private static Vector3d CalculateEccentricityVector(double gravitationalParameter, double velocityMag, double positionMag, Vector3d position, Vector3d velocity)
        {
            return (1d/gravitationalParameter) * (
                    (Pow(velocityMag,2) - gravitationalParameter/positionMag) * position
                    -
                    Vector3d.Dot(position, velocity) * velocity
                );
        }

        private static double CalculateEccentricAnomaly(double eccentricity, double trueAnomaly)
        {
            return eccentricity switch
            {
                // Elliptical:
                < 1 => 2 * Atan(Sqrt((1 - eccentricity) / (1 + eccentricity)) * Tan(trueAnomaly / 2)),
                
                // Hyperbolic:
                > 1 => 2 * Atanh(Sqrt((eccentricity - 1) / (eccentricity + 1)) * Tanh(trueAnomaly / 2)),
                
                // Circular (Parabolic):
                _ => 0
            };
        }

        private static double CalculateSemiMajorAxis(bool isParabolic, double alpha)
        {
            if (isParabolic) return double.PositiveInfinity;
            return 1d / alpha;
        }

        private static double CalculateOrbitalPeriod(bool isElliptical, double semiMajorAxis, double gravitationalParameter)
        {
            if (isElliptical)
            {
                return 2 * PI * Sqrt(
                    Pow(Abs(semiMajorAxis), 3) / gravitationalParameter
                );
            }

            return double.NaN;
        }

        // See: https://en.wikipedia.org/wiki/Stumpff_function
        private static (double, double) Stumpff(double psi)
        {
            if (psi > Accuracy)
            {
                var sqPsi = Sqrt(psi);
                var c2 = (1 - Cos(sqPsi)) / psi;
                var c3 = (sqPsi - Sin(sqPsi)) / Pow(sqPsi, 3);
                return (c2, c3);
            }

            if (psi < 0 && Abs(psi) > Accuracy)
            {
                var sqPsi = Sqrt(-psi);
                var c2 = (1 - Cosh(sqPsi)) / psi;
                var c3 = (Sinh(sqPsi) - sqPsi) / Sqrt(Pow(-psi, 3));
                return (c2, c3);
            }

            return (0.5, 1.0 / 6.0);
        }
        
        // Logic:
        private static double Cotangent(double x) => 1d / Tan(x);

        public (Vector3d, Vector3d) GetPositionAndVelocity(double deltaTime)
        {
            var p = _orbitalPeriod;
            if (!double.IsNaN(p))
            {
                deltaTime %= p;
            }

            var dotRv = Vector3d.Dot(Position, Velocity);

            var chi = 0d;
            
            if (_isElliptical)
            {
                chi = Sqrt(_gravitationalParameter) * deltaTime * _alpha;
                if (Abs(_alpha) - 1 < Accuracy)
                {
                    chi *= 0.97;
                }
            } else if (_isParabolic)
            {
                var _s = (
                    (PI / 2) - Atan(
                        deltaTime * 3 * Sqrt(
                            _gravitationalParameter / Pow(_semilatusRectum, 3)
                        )
                    )
                ) / 2.0;
                var _w = Atan(Pow(Tan(_s),1d/3d));
                chi = Sqrt(_semilatusRectum) * (2 * Cotangent(2 * _w));
            } 
            else if (_isHyperbolic)
            {
                var _num = (-2 * _gravitationalParameter * deltaTime * _alpha) / (
                    dotRv + Sign(deltaTime) * Sqrt(-_gravitationalParameter * _semiMajorAxis) *
                    (1 - _positionMag * _alpha)
                );

                if (_num < Accuracy)
                {
                    _num = Accuracy;
                }

                chi = Sign(deltaTime) * Sqrt(-_semiMajorAxis) * Log(_num);
            }
            else
            {
                return (Position, Velocity);
            }

            var iterations = 0;
            var _chi = chi + 1d;

            double c2 = 0, c3 = 0, psi = 0;

            var gpSqrt = Sqrt(_gravitationalParameter);
            
            while (Abs(chi - _chi) > Accuracy && iterations++ < 300)
            {

                psi = Pow(chi, 2) * _alpha;
                (c2, c3) = Stumpff(psi);

                var _r = Pow(chi, 2) * c2 +
                         (dotRv / gpSqrt) * chi * (1 - psi * c3) +
                         _positionMag * (1 - psi * c2);

                var chi_ = chi +
                           (
                               gpSqrt * deltaTime -
                               Pow(chi, 3) * c3 -
                               dotRv / gpSqrt * Pow(chi, 2) * c2 -
                               _positionMag * chi * (1 - psi * c3)
                           ) / _r;
                
                // Reassign for next iteration
                _chi = chi;
                chi = chi_;
            }

            if (iterations == 0)
            {
                psi = Pow(chi, 2) * _alpha;
                (c2, c3) = Stumpff(psi);
            }
            
            // At long last, we calculate the new position and velocity

            var f = 1 - (Pow(chi, 2) * c2 / _positionMag);
            var g = deltaTime - Pow(chi, 3) * c3 / Sqrt(_gravitationalParameter);

            var pos = new Vector3d(
                (f * Position.x + g * Velocity.x),
                (f * Position.y + g * Velocity.y),
                (f * Position.z + g * Velocity.z)
            );
            var _posMag = pos.magnitude;

            var g_dot = 1d - (Pow(chi, 2) * c2/_posMag);
            var f_dot = (gpSqrt * chi / (_posMag * _posMag)) * (psi * c3 - 1d);

            var vel = new Vector3d(
                f_dot * pos.x + g_dot * Velocity.x,
                f_dot * pos.y + g_dot * Velocity.y,
                f_dot * pos.z + g_dot * Velocity.z
                );

            return (pos, vel);
        }
        
        
        public (Vector3d position, Vector3d velocity) GetPositionAndVelocity2(double deltaTime)
        {
            // 1. Calculate Mean Anomaly
            var meanAnomaly = _meanMotion * deltaTime;

            // 2. Solve Kepler's Equation for Eccentric Anomaly
            double eccentricAnomaly = SolveKeplersEquation(meanAnomaly, _eccentricity);

            // 3. Calculate True Anomaly
            var trueAnomaly = CalculateTrueAnomaly(eccentricAnomaly, _eccentricity);
            
            // 4. Calculate Position and Velocity in the Perifocal Coordinate System
            //Vector3d position = CalculatePositionInPerifocal2(trueAnomaly);
            //Vector3d velocity = CalculateVelocityInPerifocal(trueAnomaly, eccentricAnomaly);
            var rDistance = _semiMajorAxis * (1 - _eccentricity * Cos(eccentricAnomaly));

            var position = rDistance * new Vector3d(
                Cos(trueAnomaly),
                Sin(trueAnomaly),
                0
            );

            var velocity = Sqrt(_gravitationalParameter * _semiMajorAxis) / rDistance * new Vector3d(
                -Sin(eccentricAnomaly),
                Sqrt(1 - Pow(_eccentricity, 2)) * Cos(eccentricAnomaly),
                0
            );
            
            
            // 5. Transform to Geocentric Equatorial Coordinate System (if required)
            // This step depends on your specific requirements
            // Assuming transformation is required and is implemented in TransformToGEC method
            position = TransformToECI(position);
            velocity = TransformToECI(velocity);

            return (position, velocity);
        }
        
        private Vector3d CalculatePositionInPerifocal(double trueAnomaly)
        {
            var r = _semiMajorAxis * (1 - _eccentricity * Cos(trueAnomaly));
            return new Vector3d(r * Cos(trueAnomaly), r * Sin(trueAnomaly), 0);
        }
        
        private Vector3d CalculatePositionInPerifocal2(double trueAnomaly)
        {
            double r = _semiMajorAxis * (1 - Pow(_eccentricity, 2)) / (1 + _eccentricity * Cos(trueAnomaly));
            return new Vector3d(r * Cos(trueAnomaly), r * Sin(trueAnomaly), 0);
        }
        
        private Vector3d CalculateVelocityInPerifocal(double trueAnomaly, double eccentricAnomaly)
        {
            var r = _semiMajorAxis * (1 - _eccentricity * Cos(eccentricAnomaly));
            var rDot = Sqrt(_gravitationalParameter / _semiMajorAxis) * _eccentricity * Sin(eccentricAnomaly);
            var trueAnomalyDot = Sqrt(_gravitationalParameter / _semiMajorAxis) * Sqrt(1 - Pow(_eccentricity, 2)) / r;

            var xDot = rDot * Cos(trueAnomaly) - r * trueAnomalyDot * Sin(trueAnomaly);
            var yDot = rDot * Sin(trueAnomaly) + r * trueAnomalyDot * Cos(trueAnomaly);

            return new Vector3d(xDot, yDot, 0);
        }

        private Vector3d TransformToECI(Vector3d vector)
        {
            // Rotation for the longitude of the ascending node
            var rotLongOfAscendingNode = QuaternionD.AngleAxis(-_longitudeOfAscendingNode * 57.2958, new Vector3d(0, 0, 1));

            // Rotation for the inclination
            var rotInclination = QuaternionD.AngleAxis(_inclination * 57.2958, new Vector3d(0, 1, 0));

            // Rotation for the argument of periapsis
            var rotArgumentOfPeriapsis = QuaternionD.AngleAxis(_argumentOfPeriapsis * 57.2958, new Vector3d(0, 0, 1));

            // Combine rotations
            var combinedRotation = rotLongOfAscendingNode * rotArgumentOfPeriapsis * rotInclination * rotLongOfAscendingNode;

            // Apply rotation to the vector
            var rotatedVector = combinedRotation * vector;

            return rotatedVector;
            return vector;
        }
        

        private Vector3d RotateAboutZAxis(Vector3d vector, double angle)
        {
            double cosAngle = Cos(angle);
            double sinAngle = Sin(angle);
            return new Vector3d(
                cosAngle * vector.x - sinAngle * vector.y,
                sinAngle * vector.x + cosAngle * vector.y,
                vector.z
            );
        }

        private Vector3d RotateAboutYAxis(Vector3d vector, double angle)
        {
            double cosAngle = Cos(angle);
            double sinAngle = Sin(angle);
            return new Vector3d(
                cosAngle * vector.x + sinAngle * vector.z,
                vector.y,
                -sinAngle * vector.x + cosAngle * vector.z
            );
        }


        
        private static double CalculateTrueAnomaly(double eccentricAnomaly, double eccentricity)
        {
            return 2 * Atan2(Sqrt(1 + eccentricity) * Sin(eccentricAnomaly / 2), Sqrt(1 - eccentricity) * Cos(eccentricAnomaly / 2));
        }
        
        private static double SolveKeplersEquation(double meanAnomaly, double eccentricity)
        {
            var result = meanAnomaly;
            var deltaE = 0.0;
            do
            {
                deltaE = (meanAnomaly + eccentricity * Sin(result) - result) / (1 - eccentricity * Cos(result));
                result += deltaE;
            }
            while (Abs(deltaE) > 1e-6); // Precision threshold

            return result;
        }
        
        private static double Kepler(double meanAnomaly, double eccentricity){
            var result = meanAnomaly;
            var breakout = 1000;
            while (breakout-- > 0)
            {
                var diff = result - eccentricity * Sin(result) - meanAnomaly;
                result -= diff / (1 - eccentricity * Cos(result));
                if (Abs(diff) <= Accuracy)
                {
                    break;
                }
            }

            return result;
        }

        public static OrbitalElements Reverse(CelestialBody parent,
            double trueAnomaly,
            double eccentricAnomaly,
            double eccentricity,
            double gravitationalParameter,
            double semiMajorAxis,
            double longOfAscendingNode,
            double inclination,
            double argumentOfPeriapsis)
        {
            var rDistance = semiMajorAxis * (1 - eccentricity * Cos(eccentricAnomaly));

            var posOrbital = rDistance * new Vector3d(
                Cos(trueAnomaly),
                Sin(trueAnomaly),
                0
            );

            var velOrbital = Sqrt(gravitationalParameter * semiMajorAxis) / rDistance * new Vector3d(
                -Sin(eccentricAnomaly),
                Sqrt(1 - Pow(eccentricity, 2)) * Cos(eccentricAnomaly),
                0
            );
            
            var vPosition = RotateOrbitalVector(posOrbital, longOfAscendingNode, inclination, argumentOfPeriapsis);
            var vVelocity = RotateOrbitalVector(velOrbital, longOfAscendingNode, inclination, argumentOfPeriapsis);
            return new OrbitalElements(parent, vPosition, vVelocity);
        }
        
        private static Vector3d RotateOrbitalVector(Vector3d vector, double longOfAscendingNode, double inclination, double argumentOfPeriapsis)
        {
            // Rotation for the longitude of the ascending node
            var rotLongOfAscendingNode = QuaternionD.AngleAxis(longOfAscendingNode * 57.2958, new Vector3d(0, 0, 1));

            // Rotation for the inclination
            var rotInclination = QuaternionD.AngleAxis(inclination * 57.2958, new Vector3d(1, 0, 0));

            // Rotation for the argument of periapsis
            var rotArgumentOfPeriapsis = QuaternionD.AngleAxis(argumentOfPeriapsis * 57.2958, new Vector3d(0, 0, 1));

            // Combine rotations
            var combinedRotation = rotArgumentOfPeriapsis * rotInclination * rotLongOfAscendingNode;

            // Apply rotation to the vector
            var rotatedVector = combinedRotation * vector;

            return rotatedVector;
        }

        
        public static OrbitalElements FromTLEs(
            double tleMeanMotion, 
            double tleMeanAnomaly, 
            double tleEccentricity, 
            double tleLongOfAscendingNode, 
            double tleInclination, 
            double tleArgumentOfPeriapsis, 
            CelestialBody parent)
        {
            var gravitationalParameter = G * parent.Mass;
            var semiMajorAxis = Cbrt(Pow((tleMeanMotion / (2 * PI)), 2) * gravitationalParameter);

            var meanAnomaly = tleMeanAnomaly * Sqrt(gravitationalParameter / Pow(semiMajorAxis, 3));
            meanAnomaly %= 2 * PI;

            var eccentricity = tleEccentricity;

            var eccentricAnomaly = Kepler(meanAnomaly, eccentricity);

            var arg1 = Sqrt(1 + eccentricity) * Sin(eccentricAnomaly / 2d);
            var arg2 = Sqrt(1 - eccentricity) * Cos(eccentricAnomaly / 2d);
            var trueAnomaly = 2 * Atan2(arg1, arg2);

            return Reverse(
                parent,
                trueAnomaly,
                eccentricAnomaly,
                eccentricity,
                gravitationalParameter,
                semiMajorAxis,
                tleLongOfAscendingNode,
                tleInclination,
                tleArgumentOfPeriapsis
            );
        }
    }
}