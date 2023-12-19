using System;
using Avionic.Math;
using Unity.Collections;
using Unity.Jobs;
using UnityEngine;

namespace Avionic.Simulation.Kepler
{
	/// <summary>
	/// Orbit data container.
	/// Also contains methods for altering and updating orbit state.
	/// </summary>
	public class OrbitalData
	{
		/// <summary>
		/// Normal of ecliptic plane.
		/// </summary>
		public static readonly Vector3d EclipticNormal = new(0, 0, 1);

		/// <summary>
		/// Up direction on ecliptic plane (y-axis on xy ecliptic plane).
		/// </summary>
		public static readonly Vector3d EclipticUp = new(0, 1, 0);

		/// <summary>
		/// Right vector on ecliptic plane (x-axis on xy ecliptic plane).
		/// </summary>
		public static readonly Vector3d EclipticRight = new(1, 0, 0);

		/// <summary>
		/// Body position relative to attractor or Focal Position.
		/// </summary>
		/// <remarks>
		/// Attractor (focus) is local center of orbit system.
		/// </remarks>
		public Vector3d Position;

		/// <summary>
		/// Magnitude of body position vector.
		/// </summary>
		public double AttractorDistance;

		/// <summary>
		/// Attractor point mass.
		/// </summary>
		public double AttractorMass;

		/// <summary>
		/// Body velocity vector relative to attractor.
		/// </summary>
		public Vector3d Velocity;

		public DateTime Epoch;

		/// <summary>
		/// Gravitational parameter of system.
		/// </summary>
		public double Mg;

		public double   SemiMinorAxis;
		public double   SemiMajorAxis;
		public double   FocalParameter;
		public double   Eccentricity;
		public double   Period;
		public double   TrueAnomaly;
		public double   MeanAnomaly;
		public double   EccentricAnomaly;
		public double   MeanMotion;
		public Vector3d Periapsis;
		public double   PeriapsisDistance;
		public Vector3d Apoapsis;
		public double   ApoapsisDistance;
		public Vector3d CenterPoint;
		public double   OrbitCompressionRatio;
		public Vector3d OrbitNormal;
		public Vector3d SemiMinorAxisBasis;
		public Vector3d SemiMajorAxisBasis;

		/// <summary>
		/// if > 0, then orbit motion is clockwise
		/// </summary>
		public double OrbitNormalDotEclipticNormal;

		public double EnergyTotal => Velocity.sqrMagnitude - 2 * Mg / AttractorDistance;

		/// <summary>
		/// The orbit inclination in radians relative to ecliptic plane.
		/// </summary>
		public double Inclination
		{
			get
			{
				var dot = Vector3d.Dot(OrbitNormal, EclipticNormal);
				return System.Math.Acos(dot);
			}
		}

		/// <summary>
		/// Ascending node longitude in radians.
		/// </summary>
		public double AscendingNodeLongitude
		{
			get
			{
				var ascNodeDir = Vector3d.Cross(EclipticNormal, OrbitNormal).normalized;
				var dot        = Vector3d.Dot(ascNodeDir, EclipticRight);
				var angle      = System.Math.Acos(dot);
				if (Vector3d.Dot(Vector3d.Cross(ascNodeDir, EclipticRight), EclipticNormal) >= 0)
				{
					angle = KeplerOrbitUtils.PI_2 - angle;
				}

				return angle;
			}
		}

		/// <summary>
		/// Angle between main orbit axis and ecliptic 0 axis in radians.
		/// </summary>
		public double ArgumentOfPerifocus
		{
			get
			{
				var ascNodeDir = Vector3d.Cross(EclipticNormal, OrbitNormal).normalized;
				var dot        = Vector3d.Dot(ascNodeDir, SemiMajorAxisBasis.normalized);
				var angle      = System.Math.Acos(dot);
				if (Vector3d.Dot(Vector3d.Cross(ascNodeDir, SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					angle = KeplerOrbitUtils.PI_2 - angle;
				}

				return angle;
			}
		}

		/// <summary>
		/// Is orbit state valid and error-free.
		/// </summary>
		/// <value>
		///   <c>true</c> if this instance is valid orbit; otherwise, <c>false</c>.
		/// </value>
		public bool IsValidOrbit =>
			Eccentricity >= 0
			&& Period > 0
			&& AttractorMass > 0;

		/// <summary>
		/// Create new orbit state without initialization. Manual orbit initialization is required.
		/// </summary>
		/// <remarks>
		/// To manually initialize orbit, fill known orbital elements and then call CalculateOrbitStateFrom... method.
		/// </remarks>
		public OrbitalData()
		{
		}

		/// <summary>
		/// Create and initialize new orbit state.
		/// </summary>
		/// <param name="position">Body local position, relative to attractor.</param>
		/// <param name="velocity">Body local velocity.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational Constant.</param>
		public OrbitalData(Vector3d position, Vector3d velocity, double attractorMass)
		{
			this.Position      = position;
			this.Velocity      = velocity;
			this.AttractorMass = attractorMass;
			CalculateOrbitStateFromOrbitalVectors();
		}

		/// <summary>
		/// Create and initialize new orbit state from orbital elements.
		/// </summary>
		/// <param name="eccentricity">Eccentricity.</param>
		/// <param name="semiMajorAxis">Main axis semi width.</param>
		/// <param name="meanAnomalyDeg">Mean anomaly in degrees.</param>
		/// <param name="inclinationDeg">Orbit inclination in degrees.</param>
		/// <param name="argOfPerifocusDeg">Orbit argument of perifocus in degrees.</param>
		/// <param name="ascendingNodeDeg">Longitude of ascending node in degrees.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational constant.</param>
		public OrbitalData(double eccentricity, double semiMajorAxis, double meanAnomalyDeg, double inclinationDeg, double argOfPerifocusDeg, double ascendingNodeDeg, double attractorMass)
		{
			this.Eccentricity  = eccentricity;
			this.SemiMajorAxis = semiMajorAxis;
			if (eccentricity < 1.0)
			{
				this.SemiMinorAxis = SemiMajorAxis * System.Math.Sqrt(1 - Eccentricity * Eccentricity);
			}
			else if (eccentricity > 1.0)
			{
				this.SemiMinorAxis = SemiMajorAxis * System.Math.Sqrt(Eccentricity * Eccentricity - 1);
			}
			else
			{
				this.SemiMajorAxis = 0;
			}

			var normal        = EclipticNormal.normalized;
			var ascendingNode = EclipticRight.normalized;

			ascendingNodeDeg %= 360;
			if (ascendingNodeDeg > 180) ascendingNodeDeg -= 360;
			inclinationDeg %= 360;
			if (inclinationDeg > 180) inclinationDeg -= 360;
			argOfPerifocusDeg %= 360;
			if (argOfPerifocusDeg > 180) argOfPerifocusDeg -= 360;

			ascendingNode = KeplerOrbitUtils.RotateVectorByAngle(ascendingNode, ascendingNodeDeg * KeplerOrbitUtils.Deg2Rad, normal).normalized;
			normal        = KeplerOrbitUtils.RotateVectorByAngle(normal,        inclinationDeg * KeplerOrbitUtils.Deg2Rad,   ascendingNode).normalized;
			Periapsis     = ascendingNode;
			Periapsis     = KeplerOrbitUtils.RotateVectorByAngle(Periapsis, argOfPerifocusDeg * KeplerOrbitUtils.Deg2Rad, normal).normalized;

			this.SemiMajorAxisBasis = Periapsis;
			this.SemiMinorAxisBasis = Vector3d.Cross(Periapsis, normal);

			this.MeanAnomaly      = meanAnomalyDeg * KeplerOrbitUtils.Deg2Rad;
			this.EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(this.MeanAnomaly, this.Eccentricity);
			this.TrueAnomaly      = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(this.EccentricAnomaly, this.Eccentricity);
			this.AttractorMass    = attractorMass;
			CalculateOrbitStateFromOrbitalElements();
		}

		/// <summary>
		/// Create and initialize new orbit state from orbital elements and main axis vectors.
		/// </summary>
		/// <param name="eccentricity">Eccentricity.</param>
		/// <param name="semiMajorAxis">Semi major axis vector.</param>
		/// <param name="semiMinorAxis">Semi minor axis vector.</param>
		/// <param name="meanAnomalyDeg">Mean anomaly in degrees.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational constant.</param>
		public OrbitalData(double eccentricity, Vector3d semiMajorAxis, Vector3d semiMinorAxis, double meanAnomalyDeg, double attractorMass)
		{
			this.Eccentricity       = eccentricity;
			this.SemiMajorAxisBasis = semiMajorAxis.normalized;
			this.SemiMinorAxisBasis = semiMinorAxis.normalized;
			this.SemiMajorAxis      = semiMajorAxis.magnitude;
			this.SemiMinorAxis      = semiMinorAxis.magnitude;

			this.MeanAnomaly      = meanAnomalyDeg * KeplerOrbitUtils.Deg2Rad;
			this.EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(this.MeanAnomaly, this.Eccentricity);
			this.TrueAnomaly      = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(this.EccentricAnomaly, this.Eccentricity);
			this.AttractorMass    = attractorMass;
			CalculateOrbitStateFromOrbitalElements();
		}

		/// <summary>
		/// Calculates full orbit state from cartesian vectors: current body position, velocity, attractor mass, and gravConstant.
		/// </summary>
		public void CalculateOrbitStateFromOrbitalVectors()
		{
			Mg                = AttractorMass * MathConstants.G;
			AttractorDistance = Position.magnitude;
			var angularMomentumVector = Vector3d.Cross(Position, Velocity);
			OrbitNormal = angularMomentumVector.normalized;
			Vector3d eccVector;
			if (OrbitNormal.sqrMagnitude < 0.99)
			{
				OrbitNormal = Vector3d.Cross(Position, EclipticUp).normalized;
				eccVector   = new Vector3d();
			}
			else
			{
				eccVector = Vector3d.Cross(Velocity, angularMomentumVector) / Mg - Position / AttractorDistance;
			}

			OrbitNormalDotEclipticNormal = Vector3d.Dot(OrbitNormal, EclipticNormal);
			FocalParameter               = angularMomentumVector.sqrMagnitude / Mg;
			Eccentricity                 = eccVector.magnitude;
			SemiMinorAxisBasis           = Vector3d.Cross(angularMomentumVector, -eccVector).normalized;
			if (SemiMinorAxisBasis.sqrMagnitude < 0.99)
			{
				SemiMinorAxisBasis = Vector3d.Cross(OrbitNormal, Position).normalized;
			}

			SemiMajorAxisBasis = Vector3d.Cross(OrbitNormal, SemiMinorAxisBasis).normalized;
			if (Eccentricity < 1.0)
			{
				OrbitCompressionRatio = 1 - Eccentricity * Eccentricity;
				SemiMajorAxis         = FocalParameter / OrbitCompressionRatio;
				SemiMinorAxis         = SemiMajorAxis * System.Math.Sqrt(OrbitCompressionRatio);
				CenterPoint           = -SemiMajorAxis * eccVector;
				var p = System.Math.Sqrt(System.Math.Pow(SemiMajorAxis, 3) / Mg);
				Period            = KeplerOrbitUtils.PI_2 * p;
				MeanMotion        = 1d / p;
				Apoapsis          = CenterPoint - SemiMajorAxisBasis * SemiMajorAxis;
				Periapsis         = CenterPoint + SemiMajorAxisBasis * SemiMajorAxis;
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance  = Apoapsis.magnitude;
				TrueAnomaly       = Vector3d.Angle(Position, SemiMajorAxisBasis) * KeplerOrbitUtils.Deg2Rad;
				if (Vector3d.Dot(Vector3d.Cross(Position, -SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}

				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(TrueAnomaly, Eccentricity);
				MeanAnomaly      = EccentricAnomaly - Eccentricity * System.Math.Sin(EccentricAnomaly);
			}
			else if (Eccentricity > 1.0)
			{
				OrbitCompressionRatio = Eccentricity * Eccentricity - 1;
				SemiMajorAxis         = FocalParameter / OrbitCompressionRatio;
				SemiMinorAxis         = SemiMajorAxis * System.Math.Sqrt(OrbitCompressionRatio);
				CenterPoint           = SemiMajorAxis * eccVector;
				Period                = double.PositiveInfinity;
				MeanMotion            = System.Math.Sqrt(Mg / System.Math.Pow(SemiMajorAxis, 3));
				Apoapsis              = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				Periapsis             = CenterPoint - SemiMajorAxisBasis * (SemiMajorAxis);
				PeriapsisDistance     = Periapsis.magnitude;
				ApoapsisDistance      = double.PositiveInfinity;
				TrueAnomaly           = Vector3d.Angle(Position, eccVector) * KeplerOrbitUtils.Deg2Rad;
				if (Vector3d.Dot(Vector3d.Cross(Position, -SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					TrueAnomaly = -TrueAnomaly;
				}

				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(TrueAnomaly, Eccentricity);
				MeanAnomaly      = System.Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
			else
			{
				OrbitCompressionRatio = 0;
				SemiMajorAxis         = 0;
				SemiMinorAxis         = 0;
				PeriapsisDistance     = angularMomentumVector.sqrMagnitude / Mg;
				CenterPoint           = new Vector3d();
				Periapsis             = -PeriapsisDistance * SemiMinorAxisBasis;
				Period                = double.PositiveInfinity;
				MeanMotion            = System.Math.Sqrt(Mg / System.Math.Pow(PeriapsisDistance, 3));
				Apoapsis              = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				ApoapsisDistance      = double.PositiveInfinity;
				TrueAnomaly           = Vector3d.Angle(Position, eccVector) * KeplerOrbitUtils.Deg2Rad;
				if (Vector3d.Dot(Vector3d.Cross(Position, -SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					TrueAnomaly = -TrueAnomaly;
				}

				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(TrueAnomaly, Eccentricity);
				MeanAnomaly      = System.Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
		}

		/// <summary>
		/// Calculates the full state of orbit from current orbital elements: eccentricity, mean anomaly, semi major and semi minor axis.
		/// </summary>
		/// <remarks>
		/// Update orbital state using known main orbital elements and basis axis vectors.
		/// Can be used for first initialization of orbit state, in this case initial data must be filled before this method call.
		/// Required initial data: eccentricity, mean anomaly, inclination, attractor mass, grav constant, all anomalies, semi minor and semi major axis vectors and magnitudes.
		/// Note that semi minor and semi major axis must be fully precalculated from inclination and argument of periapsis or another source data;
		/// </remarks>
		public void CalculateOrbitStateFromOrbitalElements()
		{
			Mg                           = AttractorMass * MathConstants.G;
			OrbitNormal                  = -Vector3d.Cross(SemiMajorAxisBasis, SemiMinorAxisBasis).normalized;
			OrbitNormalDotEclipticNormal = Vector3d.Dot(OrbitNormal, EclipticNormal);
			if (Eccentricity < 1.0)
			{
				OrbitCompressionRatio = 1 - Eccentricity * Eccentricity;
				CenterPoint           = -SemiMajorAxisBasis * SemiMajorAxis * Eccentricity;
				Period                = KeplerOrbitUtils.PI_2 * System.Math.Sqrt(System.Math.Pow(SemiMajorAxis, 3) / Mg);
				MeanMotion            = KeplerOrbitUtils.PI_2 / Period;
				Apoapsis              = CenterPoint - SemiMajorAxisBasis * SemiMajorAxis;
				Periapsis             = CenterPoint + SemiMajorAxisBasis * SemiMajorAxis;
				PeriapsisDistance     = Periapsis.magnitude;
				ApoapsisDistance      = Apoapsis.magnitude;
				// All anomalies state already preset.
			}
			else if (Eccentricity > 1.0)
			{
				CenterPoint       = SemiMajorAxisBasis * SemiMajorAxis * Eccentricity;
				Period            = double.PositiveInfinity;
				MeanMotion        = System.Math.Sqrt(Mg / System.Math.Pow(SemiMajorAxis, 3));
				Apoapsis          = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				Periapsis         = CenterPoint - SemiMajorAxisBasis * (SemiMajorAxis);
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance  = double.PositiveInfinity;
			}
			else
			{
				CenterPoint       = new Vector3d();
				Period            = double.PositiveInfinity;
				MeanMotion        = System.Math.Sqrt(Mg * 0.5 / System.Math.Pow(PeriapsisDistance, 3));
				Apoapsis          = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				PeriapsisDistance = SemiMajorAxis;
				SemiMajorAxis     = 0;
				Periapsis         = -PeriapsisDistance * SemiMajorAxisBasis;
				ApoapsisDistance  = double.PositiveInfinity;
			}

			Position = GetFocalPositionAtEccentricAnomaly(EccentricAnomaly);
			var compression = Eccentricity < 1 ? (1 - Eccentricity * Eccentricity) : (Eccentricity * Eccentricity - 1);
			FocalParameter    = SemiMajorAxis * compression;
			Velocity          = GetVelocityAtTrueAnomaly(this.TrueAnomaly);
			AttractorDistance = Position.magnitude;
		}

		/// <summary>
		/// Gets the velocity vector value at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Velocity vector.</returns>
		public Vector3d GetVelocityAtEccentricAnomaly(double eccentricAnomaly)
		{
			return GetVelocityAtTrueAnomaly(KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(eccentricAnomaly, Eccentricity));
		}

		/// <summary>
		/// Gets the velocity value at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Velocity vector.</returns>
		public Vector3d GetVelocityAtTrueAnomaly(double trueAnomaly)
		{
			if (FocalParameter <= 0)
			{
				return new Vector3d();
			}

			var sqrtMgDivP = System.Math.Sqrt(AttractorMass * MathConstants.G / FocalParameter);
			var vX         = sqrtMgDivP * (Eccentricity + System.Math.Cos(trueAnomaly));
			var vY         = sqrtMgDivP * System.Math.Sin(trueAnomaly);
			return -SemiMinorAxisBasis * vX - SemiMajorAxisBasis * vY;
		}

		/// <summary>
		/// Gets the central position at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// </remarks>
		public Vector3d GetCentralPositionAtTrueAnomaly(double trueAnomaly)
		{
			var ecc = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(trueAnomaly, Eccentricity);
			return GetCentralPositionAtEccentricAnomaly(ecc);
		}
		
		public Vector3d GetPositionAtDateTime(DateTime now) => GetPositionAtDeltaTime((now - Epoch).TotalSeconds);

		public Vector3d GetPositionAtDeltaTime(double dt)
		{
			var eccentricAnomaly = GetEccentricAnomalyForDeltaTime(dt);
			return GetCentralPositionAtEccentricAnomaly(eccentricAnomaly) + CenterPoint;
		}

		public Vector3d GetVelocityAtDateTime(DateTime now) => GetVelocityAtDeltaTime((now - Epoch).TotalSeconds);

		public Vector3d GetVelocityAtDeltaTime(double dt)
		{
			var eccentricAnomaly = GetEccentricAnomalyForDeltaTime(dt);
			return GetVelocityAtEccentricAnomaly(eccentricAnomaly);
		}
		
		/// <summary>
		/// Gets the central position at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// </remarks>
		public Vector3d GetCentralPositionAtEccentricAnomaly(double eccentricAnomaly)
		{
			if (Eccentricity < 1.0)
			{
				var result = new Vector3d(System.Math.Sin(eccentricAnomaly) * SemiMinorAxis, -System.Math.Cos(eccentricAnomaly) * SemiMajorAxis, 0);
				return -SemiMinorAxisBasis * result.x - SemiMajorAxisBasis * result.y;
			}
			
			if (Eccentricity > 1.0)
			{
				var result = new Vector3d(System.Math.Sinh(eccentricAnomaly) * SemiMinorAxis, System.Math.Cosh(eccentricAnomaly) * SemiMajorAxis, 0);
				return -SemiMinorAxisBasis * result.x - SemiMajorAxisBasis * result.y;
			}
			
			{
				var pos = new Vector3d(
					PeriapsisDistance * System.Math.Sin(eccentricAnomaly) / (1.0 + System.Math.Cos(eccentricAnomaly)),
					PeriapsisDistance * System.Math.Cos(eccentricAnomaly) / (1.0 + System.Math.Cos(eccentricAnomaly)),
					0);
				return -SemiMinorAxisBasis * pos.x + SemiMajorAxisBasis * pos.y;
			}
		}

		/// <summary>
		/// Gets the focal position at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Position relative to attractor (focus).</returns>
		public Vector3d GetFocalPositionAtEccentricAnomaly(double eccentricAnomaly)
		{
			return GetCentralPositionAtEccentricAnomaly(eccentricAnomaly) + CenterPoint;
		}

		/// <summary>
		/// Gets the focal position at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Position relative to attractor (focus).</returns>
		public Vector3d GetFocalPositionAtTrueAnomaly(double trueAnomaly)
		{
			return GetCentralPositionAtTrueAnomaly(trueAnomaly) + CenterPoint;
		}

		/// <summary>
		/// Gets the central position.
		/// </summary>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// Orbit center point is geometric center of orbit ellipse, which is the further from attractor (Focus point) the larger eccentricity is.
		/// </remarks>
		public Vector3d GetCentralPosition()
		{
			return Position - CenterPoint;
		}

		/// <summary>
		/// Gets orbit sample points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="maxDistance">The maximum distance of points.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3d[] GetOrbitPoints(int pointsCount = 50, double maxDistance = 1000d)
		{
			return GetOrbitPoints(pointsCount, new Vector3d(), maxDistance);
		}

		/// <summary>
		/// Gets orbit sample points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="origin">The origin.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3d[] GetOrbitPoints(int pointsCount, Vector3d origin, double maxDistance = 1000d)
		{
			if (pointsCount < 2)
			{
				return Array.Empty<Vector3d>();
			}

			var result = new Vector3d[pointsCount];
			if (Eccentricity < 1.0)
			{
				if (ApoapsisDistance < maxDistance)
				{
					for (var i = 0; i < pointsCount; i++)
					{
						result[i] = GetFocalPositionAtEccentricAnomaly(i * KeplerOrbitUtils.PI_2 / (pointsCount - 1d)) + origin;
					}
				}
				else if (Eccentricity > 1.0)
				{
					var maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis, PeriapsisDistance);
					for (var i = 0; i < pointsCount; i++)
					{
						result[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
				else
				{
					var maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, PeriapsisDistance, PeriapsisDistance);
					for (var i = 0; i < pointsCount; i++)
					{
						result[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
			}
			else
			{
				if (maxDistance < PeriapsisDistance)
				{
					return Array.Empty<Vector3d>();
				}

				var maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis, PeriapsisDistance);

				for (var i = 0; i < pointsCount; i++)
				{
					result[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
				}
			}

			return result;
		}

		/// <summary>
		/// Gets the orbit sample points without unnecessary memory alloc for resulting array.
		/// However, memory allocation may occur if resulting array has not correct length.
		/// </summary>
		/// <param name="orbitPoints">The orbit points.</param>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="origin">The origin.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		public void GetOrbitPointsNoAlloc(ref Vector3d[] orbitPoints, int pointsCount, Vector3d origin, double maxDistance = 1e12)
		{
			if (pointsCount < 2)
			{
				orbitPoints = Array.Empty<Vector3d>();
				return;
			}

			if (Eccentricity < 1)
			{
				if (orbitPoints == null || orbitPoints.Length != pointsCount)
				{
					orbitPoints = new Vector3d[pointsCount];
				}

				if (ApoapsisDistance < maxDistance)
				{
					for (var i = 0; i < pointsCount; i++)
					{
						orbitPoints[i] = GetFocalPositionAtEccentricAnomaly(i * KeplerOrbitUtils.PI_2 / (pointsCount - 1d)) + origin;
					}
				}
				else
				{
					var maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis, PeriapsisDistance);
					for (var i = 0; i < pointsCount; i++)
					{
						orbitPoints[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
			}
			else
			{
				if (maxDistance < PeriapsisDistance)
				{
					orbitPoints = Array.Empty<Vector3d>();
					return;
				}

				if (orbitPoints == null || orbitPoints.Length != pointsCount)
				{
					orbitPoints = new Vector3d[pointsCount];
				}

				var maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis, PeriapsisDistance);

				for (var i = 0; i < pointsCount; i++)
				{
					orbitPoints[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
				}
			}
		}

		/// <summary>
		/// Gets the ascending node of orbit.
		/// </summary>
		/// <param name="asc">The asc.</param>
		/// <returns><c>true</c> if ascending node exists, otherwise <c>false</c></returns>
		public bool GetAscendingNode(out Vector3d asc)
		{
			var    ascNodeDir  = Vector3d.Cross(OrbitNormal, EclipticNormal);
			var    s           = Vector3d.Dot(Vector3d.Cross(ascNodeDir, SemiMajorAxisBasis), OrbitNormal) >= 0;
			var    trueAnomaly = Vector3d.Angle(ascNodeDir, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
			double ecc;
			
			if (Eccentricity < 1.0)
			{
				var cosT = System.Math.Cos(trueAnomaly);
				ecc = System.Math.Acos((Eccentricity + cosT) / (1d + Eccentricity * cosT));
				if (!s)
				{
					ecc = KeplerOrbitUtils.PI_2 - ecc;
				}
			}
			else if (Eccentricity > 1.0)
			{
				trueAnomaly = Vector3d.Angle(-ascNodeDir, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
				if (trueAnomaly >= System.Math.Acos(-1d / Eccentricity))
				{
					asc = new Vector3d();
					return false;
				}

				var cosT = System.Math.Cos(trueAnomaly);
				ecc = KeplerOrbitUtils.Acosh((Eccentricity + cosT) / (1 + Eccentricity * cosT)) * (!s ? -1 : 1);
			}
			else
			{
				asc = new Vector3d();
				return false;
			}

			asc = GetFocalPositionAtEccentricAnomaly(ecc);
			return true;
		}

		/// <summary>
		/// Gets the descending node orbit.
		/// </summary>
		/// <param name="desc">The desc.</param>
		/// <returns><c>true</c> if descending node exists, otherwise <c>false</c></returns>
		public bool GetDescendingNode(out Vector3d desc)
		{
			var    norm        = Vector3d.Cross(OrbitNormal, EclipticNormal);
			var    s           = Vector3d.Dot(Vector3d.Cross(norm, SemiMajorAxisBasis), OrbitNormal) < 0;
			var    trueAnomaly = Vector3d.Angle(norm, -CenterPoint) * KeplerOrbitUtils.Deg2Rad;
			double ecc;
			
			if (Eccentricity < 1.0)
			{
				var cosT = System.Math.Cos(trueAnomaly);
				ecc = System.Math.Acos((Eccentricity + cosT) / (1d + Eccentricity * cosT));
				if (s)
				{
					ecc = KeplerOrbitUtils.PI_2 - ecc;
				}
			}
			else if (Eccentricity > 1.0)
			{
				trueAnomaly = Vector3d.Angle(norm, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
				if (trueAnomaly >= System.Math.Acos(-1d / Eccentricity))
				{
					desc = new Vector3d();
					return false;
				}

				var cosT = System.Math.Cos(trueAnomaly);
				ecc = KeplerOrbitUtils.Acosh((Eccentricity + cosT) / (1 + Eccentricity * cosT)) * (s ? -1 : 1);
			}
			else
			{
				desc = new Vector3d();
				return false;
			}

			desc = GetFocalPositionAtEccentricAnomaly(ecc);
			return true;
		}

		/// <summary>
		/// Updates the kepler orbit dynamic state (anomalies, position, velocity) by defined deltaTime.
		/// </summary>
		/// <param name="deltaTime">The delta time.</param>
		public void UpdateOrbitDataByTime(double deltaTime)
		{
			UpdateOrbitAnomaliesByTime(deltaTime);
			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}
		
		public Vector3d[] GetPositionsOverTimeBurst2(double startTime, double step, int steps)
		{
			return CalculatePositionsOverTime(startTime, step, steps);
		}
		
		public Vector3d[] GetPositionsOverTime(double startTime, double step, int steps)
		{
			var pos = new Vector3d[steps];
			for (var i = 0; i < steps; i++)
			{
				pos[i] = GetFocalPositionAtEccentricAnomaly(GetEccentricAnomalyForDeltaTime(startTime + i * steps));
			}
			return pos;
		}
		
		public void GetPositionsOverTime(Vector3d[] positions, double startTime, double step, int steps)
		{
			for (var i = 0; i < steps; i++)
			{
				positions[i] = GetFocalPositionAtEccentricAnomaly(GetEccentricAnomalyForDeltaTime(startTime + i * steps));
			}
		}

		private Vector3d[] CalculatePositionsOverTime(double startTime, double step, int steps)
		{
			Debug.Log($"E = {Eccentricity:F12}");
			var results = new NativeArray<Vector3d>(steps, Allocator.TempJob);
			if (Eccentricity <= 1e-10) // Circular
			{
				Debug.Log($"Is circular");
				var j = new PositionSolverJobCircular
				{
					StartTime = startTime,
					MeanAnomaly = MeanAnomaly,
					MeanMotion = MeanMotion,
					Step = step,
					Eccentricity = Eccentricity,
					PeriapsisDistance = PeriapsisDistance,
					SemiMinorAxisBasis = SemiMinorAxisBasis,
					SemiMajorAxisBasis = SemiMajorAxisBasis,
					Offset = CenterPoint,
					Result = results
				};
			
				j.Schedule(steps, 64).Complete();
			}
			else if (Eccentricity < 1.0) // Elliptical
			{
				var j = new PositionSolverJobElliptical
				{
					StartTime = startTime,
					MeanAnomaly = MeanAnomaly,
					MeanMotion = MeanMotion,
					Step = step,
					Eccentricity = Eccentricity,
					SemiMinorAxis = SemiMinorAxis,
					SemiMajorAxis = SemiMajorAxis,
					SemiMinorAxisBasis = SemiMinorAxisBasis,
					SemiMajorAxisBasis = SemiMajorAxisBasis,
					Offset = CenterPoint,
					Result = results
				};
				
				j.Schedule(steps, 64).Complete();
			}
			else if (Eccentricity >= 1.0) // Hyperbolic
			{
				var j = new PositionSolverJobHyperbolic
				{
					StartTime = startTime,
					MeanAnomaly = MeanAnomaly,
					MeanMotion = MeanMotion,
					Step = step,
					Eccentricity = Eccentricity,
					SemiMinorAxis = SemiMinorAxis,
					SemiMajorAxis = SemiMajorAxis,
					SemiMinorAxisBasis = SemiMinorAxisBasis,
					SemiMajorAxisBasis = SemiMajorAxisBasis,
					Offset = CenterPoint,
					Result = results
				};
				
				j.Schedule(steps, 64).Complete();
			}

			var ret = results.ToArray();
			results.Dispose();

			return ret;
		}

		public double GetEccentricAnomalyForDeltaTime(double deltaTime)
		{
			var meanAnomaly = MeanAnomaly + MeanMotion * deltaTime;
			
			if (Eccentricity >= 1.0) return KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(meanAnomaly, Eccentricity);
			
			meanAnomaly %= KeplerOrbitUtils.PI_2;
			if (meanAnomaly < 0) meanAnomaly = KeplerOrbitUtils.PI_2 - meanAnomaly;
			
			return KeplerOrbitUtils.KeplerSolver(meanAnomaly, Eccentricity);
		}

		public double GetTrueAnomalyForDeltaTime(double deltaTime)
		{
			var eccentricAnomaly = GetEccentricAnomalyForDeltaTime(deltaTime);
			
			if (Eccentricity < 1.0)
			{
				var cosE = System.Math.Cos(eccentricAnomaly);
				var trueAnomaly = System.Math.Acos((cosE - Eccentricity) / (1 - Eccentricity * cosE));
				if (trueAnomaly > System.Math.PI)
				{
					trueAnomaly = KeplerOrbitUtils.PI_2 - trueAnomaly;
				}
				return trueAnomaly;
			}
			
			if (Eccentricity > 1.0)
			{
				return System.Math.Atan2(System.Math.Sqrt(Eccentricity * Eccentricity - 1.0) * System.Math.Sinh(eccentricAnomaly), Eccentricity - System.Math.Cosh(eccentricAnomaly));
			}
			
			return eccentricAnomaly;
		}
		
		/// <summary>
		/// Updates the value of orbital anomalies by defined deltaTime.
		/// </summary>
		/// <param name="deltaTime">The delta time.</param>
		/// <remarks>
		/// Only anomalies values will be changed. 
		/// Position and velocity states needs to be updated too after this method call.
		/// </remarks>
		public void UpdateOrbitAnomaliesByTime(double deltaTime)
		{
			if (Eccentricity < 1.0)
			{
				MeanAnomaly += MeanMotion * deltaTime;
				MeanAnomaly %= KeplerOrbitUtils.PI_2;

				if (MeanAnomaly < 0)
				{
					MeanAnomaly = KeplerOrbitUtils.PI_2 - MeanAnomaly;
				}

				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				var cosE = System.Math.Cos(EccentricAnomaly);
				TrueAnomaly = System.Math.Acos((cosE - Eccentricity) / (1 - Eccentricity * cosE));
				if (MeanAnomaly > System.Math.PI)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}
			}
			else if (Eccentricity > 1.0)
			{
				MeanAnomaly      = MeanAnomaly + MeanMotion * deltaTime;
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly      = System.Math.Atan2(System.Math.Sqrt(Eccentricity * Eccentricity - 1.0) * System.Math.Sinh(EccentricAnomaly), Eccentricity - System.Math.Cosh(EccentricAnomaly));
			}
			else
			{
				MeanAnomaly      = MeanAnomaly + MeanMotion * deltaTime;
				EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(MeanAnomaly, Eccentricity);
				TrueAnomaly      = EccentricAnomaly;
			}
		}

		/// <summary>
		/// Get orbit time for elliptic orbits at current anomaly. Result is in range [0..Period].
		/// </summary>
		/// <returns>Time, corresponding to current anomaly</returns>
		public double GetCurrentOrbitTime()
		{
			if (Eccentricity < 1.0)
			{
				if (Period > 0 && Period < double.PositiveInfinity)
				{
					var anomaly = MeanAnomaly % KeplerOrbitUtils.PI_2;
					if (anomaly < 0)
					{
						anomaly = KeplerOrbitUtils.PI_2 - anomaly;
					}

					return anomaly / KeplerOrbitUtils.PI_2 * Period;
				}

				return 0.0;
			}
			else if (Eccentricity > 1.0)
			{
				var meanMotion = MeanMotion;
				if (meanMotion > 0)
				{
					return MeanAnomaly / meanMotion;
				}

				return 0.0;
			}
			else
			{
				var meanMotion = MeanMotion;
				if (meanMotion > 0)
				{
					return MeanAnomaly / meanMotion;
				}

				return 0.0;
			}
		}

		/// <summary>
		/// Updates position from anomaly state.
		/// </summary>
		public void SetPositionByCurrentAnomaly()
		{
			Position = GetFocalPositionAtEccentricAnomaly(EccentricAnomaly);
		}

		/// <summary>
		/// Sets orbit velocity, calculated by current anomaly.
		/// </summary>
		public void SetVelocityByCurrentAnomaly()
		{
			Velocity = GetVelocityAtEccentricAnomaly(EccentricAnomaly);
		}

		/// <summary>
		/// Sets the eccentricity and updates all corresponding orbit state values.
		/// </summary>
		/// <param name="e">The new eccentricity value.</param>
		public void SetEccentricity(double e)
		{
			if (!IsValidOrbit)
			{
				return;
			}

			e = System.Math.Abs(e);
			var periapsis = PeriapsisDistance; // Periapsis remains constant
			Eccentricity = e;
			var compression = Eccentricity < 1 ? (1 - Eccentricity * Eccentricity) : (Eccentricity * Eccentricity - 1);
			SemiMajorAxis  = System.Math.Abs(periapsis / (1 - Eccentricity));
			FocalParameter = SemiMajorAxis * compression;
			SemiMinorAxis  = SemiMajorAxis * System.Math.Sqrt(compression);
			CenterPoint    = SemiMajorAxis * System.Math.Abs(Eccentricity) * SemiMajorAxisBasis;
			if (Eccentricity < 1.0)
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				var cosE = System.Math.Cos(EccentricAnomaly);
				TrueAnomaly = System.Math.Acos((cosE - Eccentricity) / (1 - Eccentricity * cosE));
				if (MeanAnomaly > System.Math.PI)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}
			}
			else if (Eccentricity > 1.0)
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly      = System.Math.Atan2(System.Math.Sqrt(Eccentricity * Eccentricity - 1) * System.Math.Sinh(EccentricAnomaly), Eccentricity - System.Math.Cosh(EccentricAnomaly));
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(MeanAnomaly, Eccentricity);
				TrueAnomaly      = EccentricAnomaly;
			}

			SetVelocityByCurrentAnomaly();
			SetPositionByCurrentAnomaly();

			CalculateOrbitStateFromOrbitalVectors();
		}

		/// <summary>
		/// Sets the mean anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="m">The m.</param>
		public void SetMeanAnomaly(double m)
		{
			if (!IsValidOrbit)
			{
				return;
			}

			MeanAnomaly = m % KeplerOrbitUtils.PI_2;
			if (Eccentricity < 1.0)
			{
				if (MeanAnomaly < 0)
				{
					MeanAnomaly += KeplerOrbitUtils.PI_2;
				}

				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				TrueAnomaly      = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(EccentricAnomaly, Eccentricity);
			}
			else if (Eccentricity > 1.0)
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly      = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(EccentricAnomaly, Eccentricity);
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(MeanAnomaly, Eccentricity);
				TrueAnomaly      = EccentricAnomaly;
			}

			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Sets the true anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="t">The t.</param>
		public void SetTrueAnomaly(double t)
		{
			if (!IsValidOrbit)
			{
				return;
			}

			t %= KeplerOrbitUtils.PI_2;

			if (Eccentricity < 1.0)
			{
				if (t < 0)
				{
					t += KeplerOrbitUtils.PI_2;
				}

				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(t, Eccentricity);
				MeanAnomaly      = EccentricAnomaly - Eccentricity * System.Math.Sin(EccentricAnomaly);
			}
			else if (Eccentricity > 1.0)
			{
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(t, Eccentricity);
				MeanAnomaly      = System.Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(t, Eccentricity);
				MeanAnomaly      = KeplerOrbitUtils.ConvertEccentricToMeanAnomaly(EccentricAnomaly, Eccentricity);
			}

			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Sets the eccentric anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="e">The e.</param>
		public void SetEccentricAnomaly(double e)
		{
			if (!IsValidOrbit)
			{
				return;
			}

			e %= KeplerOrbitUtils.PI_2;

			EccentricAnomaly = e;

			if (Eccentricity < 1.0)
			{
				if (e < 0)
				{
					e = KeplerOrbitUtils.PI_2 + e;
				}

				EccentricAnomaly = e;
				TrueAnomaly      = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(e, Eccentricity);
				MeanAnomaly      = KeplerOrbitUtils.ConvertEccentricToMeanAnomaly(e, Eccentricity);
			}
			else if (Eccentricity > 1.0)
			{
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(e, Eccentricity);
				MeanAnomaly = KeplerOrbitUtils.ConvertEccentricToMeanAnomaly(e, Eccentricity);
			}
			else
			{
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(e, Eccentricity);
				MeanAnomaly = KeplerOrbitUtils.ConvertEccentricToMeanAnomaly(e, Eccentricity);
			}

			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		public object Clone()
		{
			return MemberwiseClone();
		}

		public OrbitalData CloneToUnityCoordinateSpace()
		{
			var v = -GetVelocityAtDeltaTime(0);
			var p = GetPositionAtDeltaTime(0);
			
			// Adjust these vectors to convert from ECI to Unity (Y up)
			var unityPosition = new Vector3d(p.y, p.z, p.x);
			var unityVelocity = new Vector3d(v.y, v.z, v.x);

			return new OrbitalData(unityPosition, unityVelocity, AttractorMass);
		}

		public OrbitalData CloneOrbit()
		{
			return (OrbitalData)MemberwiseClone();
		}
	}
}