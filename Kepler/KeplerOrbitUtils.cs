using System;
using Avionic.Math;

namespace Avionic.Simulation.Kepler
{
	/// <summary>
	/// Math utility methods for help in orbits calculations.
	/// </summary>
	public static class KeplerOrbitUtils
	{
		/// <summary>
		/// Two PI.
		/// </summary>
		public const double PI_2 = 6.2831853071796d;
		public const double PI = 3.14159265358979;
		public const double Deg2Rad = 0.017453292519943d;
		public const double Rad2Deg = 57.295779513082d;

		/// <summary>
		/// Regular Acosh, but without exception when out of possible range.
		/// </summary>
		/// <param name="x">The input value.</param>
		/// <returns>Calculated Acos value or 0.</returns>
		public static double Acosh(double x)
		{
			if (x < 1.0)
			{
				return 0;
			}

			return System.Math.Log(x + System.Math.Sqrt(x * x - 1.0));
		}

		/// <summary>
		/// Rotate vector around another vector (double).
		/// </summary>
		/// <param name="v">Vector to rotate.</param>
		/// <param name="angleRad">angle in radians.</param>
		/// <param name="n">normalized vector to rotate around, or normal of rotation plane.</param>
		public static Vector3d RotateVectorByAngle(Vector3d v, double angleRad, Vector3d n)
		{
			var cosT        = System.Math.Cos(angleRad);
			var sinT        = System.Math.Sin(angleRad);
			var oneMinusCos = 1f - cosT;
			// Rotation matrix:
			var a11 = oneMinusCos * n.x * n.x + cosT;
			var a12 = oneMinusCos * n.x * n.y - n.z * sinT;
			var a13 = oneMinusCos * n.x * n.z + n.y * sinT;
			var a21 = oneMinusCos * n.x * n.y + n.z * sinT;
			var a22 = oneMinusCos * n.y * n.y + cosT;
			var a23 = oneMinusCos * n.y * n.z - n.x * sinT;
			var a31 = oneMinusCos * n.x * n.z - n.y * sinT;
			var a32 = oneMinusCos * n.y * n.z + n.x * sinT;
			var a33 = oneMinusCos * n.z * n.z + cosT;
			return new Vector3d(
				v.x * a11 + v.y * a12 + v.z * a13,
				v.x * a21 + v.y * a22 + v.z * a23,
				v.x * a31 + v.y * a32 + v.z * a33
			);
		}

		public static double Abs(double a)
		{
			return a < 0 ? -a : a;
		}

		public static float Abs(float a)
		{
			return a < 0 ? -a : a;
		}

		/// <summary>
		/// Calculate velocity vector for circle orbit.
		/// </summary>
		public static Vector3d CalcCircleOrbitVelocity(Vector3d attractorPos, Vector3d bodyPos, double attractorMass, Vector3d orbitNormal, double gConst)
		{
			var distanceVector = bodyPos - attractorPos;
			var dist           = distanceVector.magnitude;
			var mg             = attractorMass * gConst;
			var vScalar        = System.Math.Sqrt(mg / dist);
			return Vector3d.Cross(distanceVector, -orbitNormal).normalized * vScalar;
		}

		/// <summary>
		/// Calculates the center of mass.
		/// </summary>
		/// <param name="pos1">The posistion 1.</param>
		/// <param name="mass1">The mass 1.</param>
		/// <param name="pos2">The position 2.</param>
		/// <param name="mass2">The mass 2.</param>
		/// <returns>Center of mass postion vector.</returns>
		public static Vector3d CalcCenterOfMass(Vector3d pos1, double mass1, Vector3d pos2, double mass2)
		{
			return ((pos1 * mass1) + (pos2 * mass2)) / (mass1 + mass2);
		}

		/// <summary>
		/// Converts the eccentric to true anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>True anomaly in radians.</returns>
		public static double ConvertEccentricToTrueAnomaly(double eccentricAnomaly, double eccentricity)
		{
			if (eccentricity < 1.0)
			{
				var cosE  = System.Math.Cos(eccentricAnomaly);
				var tAnom = System.Math.Acos((cosE - eccentricity) / (1d - eccentricity * cosE));
				if (eccentricAnomaly > PI)
				{
					tAnom = PI_2 - tAnom;
				}

				return tAnom;
			}
			else if (eccentricity > 1.0)
			{
				var tAnom = System.Math.Atan2(
					System.Math.Sqrt(eccentricity * eccentricity - 1d) * System.Math.Sinh(eccentricAnomaly),
					eccentricity - System.Math.Cosh(eccentricAnomaly)
				);
				return tAnom;
			}
			else
			{
				return eccentricAnomaly;
			}
		}

		/// <summary>
		/// Converts the true to eccentric anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>Eccentric anomaly in radians.</returns>
		public static double ConvertTrueToEccentricAnomaly(double trueAnomaly, double eccentricity)
		{
			if (double.IsNaN(eccentricity) || double.IsInfinity(eccentricity))
			{
				return trueAnomaly;
			}

			trueAnomaly %= PI_2;
			if (eccentricity < 1.0)
			{
				if (trueAnomaly < 0)
				{
					trueAnomaly += PI_2;
				}

				var cosT2   = System.Math.Cos(trueAnomaly);
				var eccAnom = System.Math.Acos((eccentricity + cosT2) / (1d + eccentricity * cosT2));
				if (trueAnomaly > System.Math.PI)
				{
					eccAnom = PI_2 - eccAnom;
				}

				return eccAnom;
			}
			else if (eccentricity > 1.0)
			{
				var cosT    = System.Math.Cos(trueAnomaly);
				var eccAnom = Acosh((eccentricity + cosT) / (1d + eccentricity * cosT)) * System.Math.Sign(trueAnomaly);
				return eccAnom;
			}
			else
			{
				// For parabolic trajectories
				// there is no Eccentric anomaly defined,
				// because 'True anomaly' to 'Time' relation can be resolved analytically.
				return trueAnomaly;
			}
		}

		/// <summary>
		/// Converts the mean to eccentric anomaly.
		/// </summary>
		/// <param name="meanAnomaly">The mean anomaly.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>Eccentric anomaly in radians.</returns>
		public static double ConvertMeanToEccentricAnomaly(double meanAnomaly, double eccentricity)
		{
			if (eccentricity < 1.0)
			{
				return KeplerSolver(meanAnomaly, eccentricity);
			}
			
			if (eccentricity > 1.0)
			{
				return KeplerSolverHyperbolicCase(meanAnomaly, eccentricity);
			}
			
			{
				var m   = meanAnomaly * 2;
				var v   = 12d * m + 4d * System.Math.Sqrt(4d + 9d * m * m);
				var pow = System.Math.Pow(v, 1d / 3d);
				var t   = 0.5 * pow - 2 / pow;
				return 2 * System.Math.Atan(t);
			}
		}

		/// <summary>
		/// Converts the eccentric to mean anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>Mean anomaly in radians.</returns>
		public static double ConvertEccentricToMeanAnomaly(double eccentricAnomaly, double eccentricity)
		{
			if (eccentricity < 1.0)
			{
				return eccentricAnomaly - eccentricity * System.Math.Sin(eccentricAnomaly);
			}
			else if (eccentricity > 1.0)
			{
				return System.Math.Sinh(eccentricAnomaly) * eccentricity - eccentricAnomaly;
			}
			else
			{
				var t = System.Math.Tan(eccentricAnomaly * 0.5);
				return (t + t * t * t / 3d) * 0.5d;
			}
		}

		/// <summary>
		/// Gets the True anomaly value from current distance from the focus (attractor).
		/// </summary>
		/// <param name="distance">The distance from attractor.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <param name="semiMajorAxis">The semi major axis.</param>
		/// <param name="periapsisDistance">The periapsis distance value.</param>
		/// <returns>True anomaly in radians.</returns>
		public static double CalcTrueAnomalyForDistance(double distance, double eccentricity, double semiMajorAxis, double periapsisDistance)
		{
			if (eccentricity < 1.0)
			{
				return System.Math.Acos((semiMajorAxis * (1d - eccentricity * eccentricity) - distance) / (distance * eccentricity));
			}
			else if (eccentricity > 1.0)
			{
				return System.Math.Acos((semiMajorAxis * (eccentricity * eccentricity - 1d) - distance) / (distance * eccentricity));
			}
			else
			{
				return System.Math.Acos((periapsisDistance / distance) - 1d);
			}
		}

		/// <summary>
		/// A classic Kepler solver.
		/// </summary>
		/// <param name="meanAnomaly">The mean anomaly in radians.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>Eccentric anomaly in radians.</returns>
		/// <remarks>
		/// One stable method.
		/// </remarks>
		public static double KeplerSolver(double meanAnomaly, double eccentricity)
		{
			// Iterations count range from 2 to 6 when eccentricity is in range from 0 to 1.
			var iterations = (int)(System.Math.Ceiling((eccentricity + 0.7d) * 1.25d)) << 1;
			var m = meanAnomaly;
			double esinE;
			double ecosE;
			double deltaE;
			double n;
			for (var i = 0; i < iterations; i++)
			{
				esinE  =  eccentricity * System.Math.Sin(m);
				ecosE  =  eccentricity * System.Math.Cos(m);
				deltaE =  m - esinE - meanAnomaly;
				n      =  1.0 - ecosE;
				m      += -5d * deltaE / (n + System.Math.Sign(n) * System.Math.Sqrt(System.Math.Abs(16d * n * n - 20d * deltaE * esinE)));
			}

			return m;
		}

		/// <summary>
		/// Kepler solver for hyperbolic case.
		/// </summary>
		/// <param name="meanAnomaly">The mean anomaly.</param>
		/// <param name="eccentricity">The eccentricity.</param>
		/// <returns>Eccentric anomaly in radians.</returns>
		public static double KeplerSolverHyperbolicCase(double meanAnomaly, double eccentricity)
		{
			var delta = 1d;

			// Adjusting Danby's guess for negative mean anomalies.
			var signM = meanAnomaly >= 0 ? 1d : -1d;
			var f = System.Math.Log(2d * System.Math.Abs(meanAnomaly) / eccentricity + 1.8d) * signM;
			if (double.IsNaN(f) || double.IsInfinity(f))
			{
				return meanAnomaly;
			}

			while (System.Math.Abs(delta) > 1e-8)
			{
				delta = (eccentricity * System.Math.Sinh(f) - f - meanAnomaly) / (eccentricity * System.Math.Cosh(f) - 1d);
				f -= delta;
			}

			return f;
		}
	}
}