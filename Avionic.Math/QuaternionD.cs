using System;
using System.Runtime.Serialization;
using System.Xml.Serialization;
using UnityEngine;
using UnityEngine.Internal;
// ReSharper disable InconsistentNaming

namespace Avionic.Math
{
	/// <summary>
	/// Quaternions are used to represent rotations.
	/// A custom completely managed implementation of UnityEngine.Quaternion
	/// Base is decompiled UnityEngine.Quaternion
	/// Doesn't implement methods marked Obsolete
	/// Does implicit coversions to and from UnityEngine.Quaternion
	///
	/// Uses code from:
	/// https://raw.githubusercontent.com/mono/opentk/master/Source/OpenTK/Math/Quaternion.cs
	/// http://answers.unity3d.com/questions/467614/what-is-the-source-code-of-quaternionlookrotation.html
	/// http://stackoverflow.com/questions/12088610/conversion-between-euler-quaternion-like-in-unity3d-engine
	/// http://stackoverflow.com/questions/11492299/quaternion-to-euler-angles-algorithm-how-to-convert-to-y-up-and-between-ha
	///
	/// Version: lpmitchell 2022-10 (author yyyy-MM)
	/// License: ODC Public Domain Dedication & License 1.0 (PDDL-1.0) https://tldrlegal.com/license/odc-public-domain-dedication-&-license-1.0-(pddl-1.0)
	/// </summary>
	[Serializable]
	[DataContract]
	public struct QuaternionD : IEquatable<QuaternionD>
	{
		public const double RadToDeg = 180.0 / System.Math.PI;
		public const double DegToRad = System.Math.PI / 180.0;

		public const double kEpsilon = 1E-06d; // should probably be used in the 0 tests in LookRotation or Slerp

		[XmlIgnore]
		public Vector3d xyz
		{
			set
			{
				x = value.x;
				y = value.y;
				z = value.z;
			}
			get
			{
				return new Vector3d(x, y, z);
			}
		}
		/// <summary>
		///   <para>X component of the Quaternion. Don't modify this directly unless you know quaternions inside out.</para>
		/// </summary>
		[DataMember(Order = 1)]
		public double x;
		/// <summary>
		///   <para>Y component of the Quaternion. Don't modify this directly unless you know quaternions inside out.</para>
		/// </summary>
		[DataMember(Order = 2)]
		public double y;
		/// <summary>
		///   <para>Z component of the Quaternion. Don't modify this directly unless you know quaternions inside out.</para>
		/// </summary>
		[DataMember(Order = 3)]
		public double z;
		/// <summary>
		///   <para>W component of the Quaternion. Don't modify this directly unless you know quaternions inside out.</para>
		/// </summary>
		[DataMember(Order = 4)]
		public double w;

		[XmlIgnore]
		public double this[int index]
		{
			get
			{
				switch (index)
				{
					case 0:
						return x;
					case 1:
						return y;
					case 2:
						return z;
					case 3:
						return w;
					default:
						throw new IndexOutOfRangeException("Invalid Quaternion index: " + index + ", can use only 0,1,2,3");
				}
			}
			set
			{
				switch (index)
				{
					case 0:
						x = value;
						break;
					case 1:
						y = value;
						break;
					case 2:
						z = value;
						break;
					case 3:
						w = value;
						break;
					default:
						throw new IndexOutOfRangeException("Invalid Quaternion index: " + index + ", can use only 0,1,2,3");
				}
			}
		}
		/// <summary>
		///   <para>The identity rotation (RO).</para>
		/// </summary>
		[XmlIgnore]
		public static QuaternionD identity
		{
			get
			{
				return new QuaternionD(0, 0, 0, 1);
			}
		}
		/// <summary>
		///   <para>Returns the euler angle representation of the rotation.</para>
		/// </summary>
		[XmlIgnore]
		public Vector3d eulerAngles
		{
			get
			{
				return ToEulerRad(this) * RadToDeg;
			}
			set
			{
				this = FromEulerRad(value * DegToRad);
			}
		}
		/// <summary>
		/// Gets the length (magnitude) of the quaternion.
		/// </summary>
		/// <seealso cref="LengthSquared"/>
		[XmlIgnore]
		public double Length
		{
			get
			{
				return System.Math.Sqrt(x * x + y * y + z * z + w * w);
			}
		}

		/// <summary>
		/// Gets the square of the quaternion length (magnitude).
		/// </summary>
		[XmlIgnore]
		public double LengthSquared
		{
			get
			{
				return x * x + y * y + z * z + w * w;
			}
		}
		/// <summary>
		///   <para>Constructs new MyQuaternion with given x,y,z,w components.</para>
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <param name="z"></param>
		/// <param name="w"></param>
		public QuaternionD(double x, double y, double z, double w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}
		/// <summary>
		/// Construct a new MyQuaternion from vector and w components
		/// </summary>
		/// <param name="v">The vector part</param>
		/// <param name="w">The w part</param>
		public QuaternionD(Vector3d v, double w)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			this.w = w;
		}
		/// <summary>
		///   <para>Set x, y, z and w components of an existing MyQuaternion.</para>
		/// </summary>
		/// <param name="new_x"></param>
		/// <param name="new_y"></param>
		/// <param name="new_z"></param>
		/// <param name="new_w"></param>
		public void Set(double new_x, double new_y, double new_z, double new_w)
		{
			x = new_x;
			y = new_y;
			z = new_z;
			w = new_w;
		}
		/// <summary>
		/// Scales the MyQuaternion to unit length.
		/// </summary>
		public void Normalize()
		{
			var scale = 1.0 / Length;
			xyz *= scale;
			w *= scale;
		}
		/// <summary>
		/// Scale the given quaternion to unit length
		/// </summary>
		/// <param name="q">The quaternion to normalize</param>
		/// <returns>The normalized quaternion</returns>
		public static QuaternionD Normalize(QuaternionD q)
		{
			QuaternionD result;
			Normalize(ref q, out result);
			return result;
		}
		/// <summary>
		/// Scale the given quaternion to unit length
		/// </summary>
		/// <param name="q">The quaternion to normalize</param>
		/// <param name="result">The normalized quaternion</param>
		public static void Normalize(ref QuaternionD q, out QuaternionD result)
		{
			var scale = 1.0 / q.Length;
			result = new QuaternionD(q.xyz * scale, q.w * scale);
		}
		/// <summary>
		///   <para>The dot product between two rotations.</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public static double Dot(QuaternionD a, QuaternionD b)
		{
			return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
		}
		/// <summary>
		///   <para>Creates a rotation which rotates /angle/ degrees around /axis/.</para>
		/// </summary>
		/// <param name="angle"></param>
		/// <param name="axis"></param>
		public static QuaternionD AngleAxis(double angle, Vector3d axis)
		{
			return AngleAxis(angle, ref axis);
		}
		private static QuaternionD AngleAxis(double degress, ref Vector3d axis)
		{
			if (axis.sqrMagnitude == 0.0)
				return identity;

			var result = identity;
			var radians = degress * DegToRad;
			radians *= 0.5;
			axis = axis.normalized;
			axis = axis * System.Math.Sin(radians);
			result.x = axis.x;
			result.y = axis.y;
			result.z = axis.z;
			result.w = System.Math.Cos(radians);

			return Normalize(result);
		}
		public void ToAngleAxis(out double angle, out Vector3d axis)
		{
			ToAxisAngleRad(this, out axis, out angle);
			angle *= RadToDeg;
		}
		/// <summary>
		///   <para>Creates a rotation which rotates from /fromDirection/ to /toDirection/.</para>
		/// </summary>
		/// <param name="fromDirection"></param>
		/// <param name="toDirection"></param>
		public static QuaternionD FromToRotation(Vector3d fromDirection, Vector3d toDirection)
		{
			return RotateTowards(LookRotation(fromDirection), LookRotation(toDirection), double.MaxValue);
		}
		/// <summary>
		///   <para>Creates a rotation which rotates from /fromDirection/ to /toDirection/.</para>
		/// </summary>
		/// <param name="fromDirection"></param>
		/// <param name="toDirection"></param>
		public void SetFromToRotation(Vector3d fromDirection, Vector3d toDirection)
		{
			this = FromToRotation(fromDirection, toDirection);
		}
		/// <summary>
		///   <para>Creates a rotation with the specified /forward/ and /upwards/ directions.</para>
		/// </summary>
		/// <param name="forward">The direction to look in.</param>
		/// <param name="upwards">The vector that defines in which direction up is.</param>
		public static QuaternionD LookRotation(Vector3d forward, [DefaultValue("Vector3d.up")] Vector3d upwards)
		{
			return LookRotation(ref forward, ref upwards);
		}
		public static QuaternionD LookRotation(Vector3d forward)
		{
			var up = new Vector3d(0,1,0);
			return LookRotation(ref forward, ref up);
		}
		// from http://answers.unity3d.com/questions/467614/what-is-the-source-code-of-quaternionlookrotation.html
		private static QuaternionD LookRotation(ref Vector3d forward, ref Vector3d up)
		{
			
			forward = forward.normalized;
			var right = Vector3d.Cross(up, forward).normalized;
			up = Vector3d.Cross(forward, right);
			var m00 = right.x;
			var m01 = right.y;
			var m02 = right.z;
			var m10 = up.x;
			var m11 = up.y;
			var m12 = up.z;
			var m20 = forward.x;
			var m21 = forward.y;
			var m22 = forward.z;


			var num8 = (m00 + m11) + m22;
			var quaternion = new QuaternionD();
			if (num8 > 0)
			{
				var num = System.Math.Sqrt(num8 + 1f);
				quaternion.w = num * 0.5;
				num = 0.5 / num;
				quaternion.x = (m12 - m21) * num;
				quaternion.y = (m20 - m02) * num;
				quaternion.z = (m01 - m10) * num;
				return quaternion;
			}
			if ((m00 >= m11) && (m00 >= m22))
			{
				var num7 = System.Math.Sqrt(((1f + m00) - m11) - m22);
				var num4 = 0.5 / num7;
				quaternion.x = 0.5 * num7;
				quaternion.y = (m01 + m10) * num4;
				quaternion.z = (m02 + m20) * num4;
				quaternion.w = (m12 - m21) * num4;
				return quaternion;
			}
			if (m11 > m22)
			{
				var num6 = System.Math.Sqrt(((1f + m11) - m00) - m22);
				var num3 = 0.5 / num6;
				quaternion.x = (m10 + m01) * num3;
				quaternion.y = 0.5 * num6;
				quaternion.z = (m21 + m12) * num3;
				quaternion.w = (m20 - m02) * num3;
				return quaternion;
			}
			var num5 = System.Math.Sqrt(((1f + m22) - m00) - m11);
			var num2 = 0.5 / num5;
			quaternion.x = (m20 + m02) * num2;
			quaternion.y = (m21 + m12) * num2;
			quaternion.z = 0.5 * num5;
			quaternion.w = (m01 - m10) * num2;
			return quaternion;
		}
		public void SetLookRotation(Vector3d view)
		{
			var up = new Vector3d(0,1,0);
			SetLookRotation(view, up);
		}
		/// <summary>
		///   <para>Creates a rotation with the specified /forward/ and /upwards/ directions.</para>
		/// </summary>
		/// <param name="view">The direction to look in.</param>
		/// <param name="up">The vector that defines in which direction up is.</param>
		public void SetLookRotation(Vector3d view, [DefaultValue("Vector3d.up")] Vector3d up)
		{
			this = LookRotation(view, up);
		}
		/// <summary>
		///   <para>Spherically interpolates between /a/ and /b/ by t. The parameter /t/ is clamped to the range [0, 1].</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <param name="t"></param>
		public static QuaternionD Slerp(QuaternionD a, QuaternionD b, double t)
		{
			return Slerp(ref a, ref b, t);
		}
		private static QuaternionD Slerp(ref QuaternionD a, ref QuaternionD b, double t)
		{
			if (t > 1) t = 1;
			if (t < 0) t = 0;
			return SlerpUnclamped(ref a, ref b, t);
		}
		/// <summary>
		///   <para>Spherically interpolates between /a/ and /b/ by t. The parameter /t/ is not clamped.</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <param name="t"></param>
		public static QuaternionD SlerpUnclamped(QuaternionD a, QuaternionD b, double t)
		{
			return SlerpUnclamped(ref a, ref b, t);
		}
		private static QuaternionD SlerpUnclamped(ref QuaternionD a, ref QuaternionD b, double t)
		{
			// if either input is zero, return the other.
			if (a.LengthSquared == 0.0)
			{
				if (b.LengthSquared == 0.0)
				{
					return identity;
				}
				return b;
			}
			else if (b.LengthSquared == 0.0)
			{
				return a;
			}
			
			var cosHalfAngle = a.w * b.w + Vector3d.Dot(a.xyz, b.xyz);

			if (cosHalfAngle >= 1.0 || cosHalfAngle <= -1.0)
			{
				// angle = 0.0, so just return one input.
				return a;
			}
			else if (cosHalfAngle < 0.0)
			{
				b.xyz = -b.xyz;
				b.w = -b.w;
				cosHalfAngle = -cosHalfAngle;
			}

			double blendA;
			double blendB;
			if (cosHalfAngle < 0.99f)
			{
				// do proper slerp for big angles
				var halfAngle = System.Math.Acos(cosHalfAngle);
				var sinHalfAngle = System.Math.Sin(halfAngle);
				var oneOverSinHalfAngle = 1.0 / sinHalfAngle;
				blendA = System.Math.Sin(halfAngle * (1.0 - t)) * oneOverSinHalfAngle;
				blendB = System.Math.Sin(halfAngle * t) * oneOverSinHalfAngle;
			}
			else
			{
				// do lerp if angle is really small.
				blendA = 1.0 - t;
				blendB = t;
			}

			var result = new QuaternionD(blendA * a.xyz + blendB * b.xyz, blendA * a.w + blendB * b.w);
			if (result.LengthSquared > 0.0)
				return Normalize(result);
			else
				return identity;
		}
		/// <summary>
		///   <para>Interpolates between /a/ and /b/ by /t/ and normalizes the result afterwards. The parameter /t/ is clamped to the range [0, 1].</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <param name="t"></param>
		public static QuaternionD Lerp(QuaternionD a, QuaternionD b, double t)
		{
			if (t > 1) t = 1;
			if (t < 0) t = 0;
			return Slerp(ref a, ref b, t); // TODO: use lerp not slerp, "Because quaternion works in 4D. Rotation in 4D are linear" ???
		}
		/// <summary>
		///   <para>Interpolates between /a/ and /b/ by /t/ and normalizes the result afterwards. The parameter /t/ is not clamped.</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <param name="t"></param>
		public static QuaternionD LerpUnclamped(QuaternionD a, QuaternionD b, double t)
		{
			return Slerp(ref a, ref b, t);
		}
		/// <summary>
		///   <para>Rotates a rotation /from/ towards /to/.</para>
		/// </summary>
		/// <param name="from"></param>
		/// <param name="to"></param>
		/// <param name="maxDegreesDelta"></param>
		public static QuaternionD RotateTowards(QuaternionD from, QuaternionD to, double maxDegreesDelta)
		{
			var num = Angle(from, to);
			if (num == 0)
			{
				return to;
			}
			var t = System.Math.Min(1f, maxDegreesDelta / num);
			return SlerpUnclamped(from, to, t);
		}
		/// <summary>
		///   <para>Returns the Inverse of /rotation/.</para>
		/// </summary>
		/// <param name="rotation"></param>
		public static QuaternionD Inverse(QuaternionD rotation)
		{
			var lengthSq = rotation.LengthSquared;
			if (lengthSq == 0.0) return rotation;
			
			var i = 1.0 / lengthSq;
			return new QuaternionD(rotation.xyz * -i, rotation.w * i);
		}
		/// <summary>
		///   <para>Returns a nicely formatted string of the MyQuaternion.</para>
		/// </summary>
		/// <param name="format"></param>
		public override string ToString()
		{
			return string.Format("({0:F1}, {1:F1}, {2:F1}, {3:F1})", x, y, z, w);
		}
		/// <summary>
		///   <para>Returns a nicely formatted string of the MyQuaternion.</para>
		/// </summary>
		/// <param name="format"></param>
		public string ToString(string format)
		{
			return string.Format("({0}, {1}, {2}, {3})", x.ToString(format), y.ToString(format), z.ToString(format), w.ToString(format));
		}
		/// <summary>
		///   <para>Returns the angle in degrees between two rotations /a/ and /b/.</para>
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public static double Angle(QuaternionD a, QuaternionD b)
		{
			var f = Dot(a, b);
			return System.Math.Acos(System.Math.Min(System.Math.Abs(f), 1)) * 2 * RadToDeg;
		}
		/// <summary>
		///   <para>Returns a rotation that rotates z degrees around the z axis, x degrees around the x axis, and y degrees around the y axis (in that order).</para>
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <param name="z"></param>
		public static QuaternionD Euler(double x, double y, double z)
		{
			return FromEulerRad(new Vector3d(x, y, z) * DegToRad);
		}
		/// <summary>
		///   <para>Returns a rotation that rotates z degrees around the z axis, x degrees around the x axis, and y degrees around the y axis (in that order).</para>
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <param name="z"></param>
		public static QuaternionD EulerRad(double x, double y, double z)
		{
			return FromEulerRad(new Vector3d(x, y, z));
		}
		/// <summary>
		///   <para>Returns a rotation that rotates z degrees around the z axis, x degrees around the x axis, and y degrees around the y axis (in that order).</para>
		/// </summary>
		/// <param name="euler"></param>
		public static QuaternionD Euler(Vector3d euler)
		{
			return FromEulerRad(euler * DegToRad);
		}
		// from http://stackoverflow.com/questions/12088610/conversion-between-euler-quaternion-like-in-unity3d-engine
		private static Vector3d ToEulerRad(QuaternionD rotation)
		{
			var sqw = rotation.w * rotation.w;
			var sqx = rotation.x * rotation.x;
			var sqy = rotation.y * rotation.y;
			var sqz = rotation.z * rotation.z;
			var unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
			var test = rotation.x * rotation.w - rotation.y * rotation.z;
			Vector3d v;

			if (test > 0.4995 * unit)
			{ // singularity at north pole
				v.y = 2 * System.Math.Atan2(rotation.y, rotation.x);
				v.x = System.Math.PI / 2;
				v.z = 0;
				return NormalizeAngles(v * RadToDeg);
			}
			if (test < -0.4995 * unit)
			{ // singularity at south pole
				v.y = -2 * System.Math.Atan2(rotation.y, rotation.x);
				v.x = -System.Math.PI / 2;
				v.z = 0;
				return NormalizeAngles(v * RadToDeg);
			}
			var q = new QuaternionD(rotation.w, rotation.z, rotation.x, rotation.y);
			v.y = System.Math.Atan2(2 * q.x * q.w + 2 * q.y * q.z, 1 - 2 * (q.z * q.z + q.w * q.w));     // Yaw
			v.x = System.Math.Asin(2 * (q.x * q.z - q.w * q.y));                             // Pitch
			v.z = System.Math.Atan2(2 * q.x * q.y + 2 * q.z * q.w, 1 - 2 * (q.y * q.y + q.z * q.z));      // Roll
			return NormalizeAngles(v * RadToDeg);
		}
		private static Vector3d NormalizeAngles(Vector3d angles)
		{
			angles.x = NormalizeAngle(angles.x);
			angles.y = NormalizeAngle(angles.y);
			angles.z = NormalizeAngle(angles.z);
			return angles;
		}
		private static double NormalizeAngle(double angle)
		{
			while (angle > 360)
				angle -= 360;
			while (angle < 0)
				angle += 360;
			return angle;
		}
		// from http://stackoverflow.com/questions/11492299/quaternion-to-euler-angles-algorithm-how-to-convert-to-y-up-and-between-ha
		private static QuaternionD FromEulerRad(Vector3d euler)
		{
			var yaw = euler.x;
			var pitch = euler.y;
			var roll = euler.z;
			var rollOver2 = roll * 0.5;
			var sinRollOver2 = System.Math.Sin(rollOver2);
			var cosRollOver2 = System.Math.Cos(rollOver2);
			var pitchOver2 = pitch * 0.5;
			var sinPitchOver2 = System.Math.Sin(pitchOver2);
			var cosPitchOver2 = System.Math.Cos(pitchOver2);
			var yawOver2 = yaw * 0.5;
			var sinYawOver2 = System.Math.Sin(yawOver2);
			var cosYawOver2 = System.Math.Cos(yawOver2);
			QuaternionD result;
			result.x = cosYawOver2 * cosPitchOver2 * cosRollOver2 + sinYawOver2 * sinPitchOver2 * sinRollOver2;
			result.y = cosYawOver2 * cosPitchOver2 * sinRollOver2 - sinYawOver2 * sinPitchOver2 * cosRollOver2;
			result.z = cosYawOver2 * sinPitchOver2 * cosRollOver2 + sinYawOver2 * cosPitchOver2 * sinRollOver2;
			result.w = sinYawOver2 * cosPitchOver2 * cosRollOver2 - cosYawOver2 * sinPitchOver2 * sinRollOver2;
			return result;

		}
		private static void ToAxisAngleRad(QuaternionD q, out Vector3d axis, out double angle)
		{
			if (System.Math.Abs(q.w) > 1.0)
				q.Normalize();
			angle = 2.0 * System.Math.Acos(q.w); // angle
			var den = System.Math.Sqrt(1.0 - q.w * q.w);
			if (den > 0.0001)
			{
				axis = q.xyz / den;
			}
			else
			{
				// This occurs when the angle is zero. 
				// Not a problem: just set an arbitrary normalized axis.
				axis = new Vector3d(1, 0, 0);
			}
		}
		#region Obsolete methods
		/*
	[Obsolete("Use MyQuaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public static MyQuaternion EulerRotation(double x, double y, double z)
	{
		return MyQuaternion.Internal_FromEulerRad(new Vector3d(x, y, z));
	}
	[Obsolete("Use MyQuaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public static MyQuaternion EulerRotation(Vector3d euler)
	{
		return MyQuaternion.Internal_FromEulerRad(euler);
	}
	[Obsolete("Use MyQuaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public void SetEulerRotation(double x, double y, double z)
	{
		this = Quaternion.Internal_FromEulerRad(new Vector3d(x, y, z));
	}
	[Obsolete("Use Quaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public void SetEulerRotation(Vector3d euler)
	{
		this = Quaternion.Internal_FromEulerRad(euler);
	}
	[Obsolete("Use Quaternion.eulerAngles instead. This function was deprecated because it uses radians instead of degrees")]
	public Vector3d ToEuler()
	{
		return Quaternion.Internal_ToEulerRad(this);
	}
	[Obsolete("Use Quaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public static Quaternion EulerAngles(double x, double y, double z)
	{
		return Quaternion.Internal_FromEulerRad(new Vector3d(x, y, z));
	}
	[Obsolete("Use Quaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public static Quaternion EulerAngles(Vector3d euler)
	{
		return Quaternion.Internal_FromEulerRad(euler);
	}
	[Obsolete("Use Quaternion.ToAngleAxis instead. This function was deprecated because it uses radians instead of degrees")]
	public void ToAxisAngle(out Vector3d axis, out double angle)
	{
		Quaternion.Internal_ToAxisAngleRad(this, out axis, out angle);
	}
	[Obsolete("Use Quaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public void SetEulerAngles(double x, double y, double z)
	{
		this.SetEulerRotation(new Vector3d(x, y, z));
	}
	[Obsolete("Use Quaternion.Euler instead. This function was deprecated because it uses radians instead of degrees")]
	public void SetEulerAngles(Vector3d euler)
	{
		this = Quaternion.EulerRotation(euler);
	}
	[Obsolete("Use Quaternion.eulerAngles instead. This function was deprecated because it uses radians instead of degrees")]
	public static Vector3d ToEulerAngles(Quaternion rotation)
	{
		return Quaternion.Internal_ToEulerRad(rotation);
	}
	[Obsolete("Use Quaternion.eulerAngles instead. This function was deprecated because it uses radians instead of degrees")]
	public Vector3d ToEulerAngles()
	{
		return Quaternion.Internal_ToEulerRad(this);
	}
	[Obsolete("Use Quaternion.AngleAxis instead. This function was deprecated because it uses radians instead of degrees")]
	public static Quaternion AxisAngle(Vector3d axis, double angle)
	{
		return Quaternion.INTERNAL_CALL_AxisAngle(ref axis, angle);
	}

	private static Quaternion INTERNAL_CALL_AxisAngle(ref Vector3d axis, double angle)
	{

	}
	[Obsolete("Use Quaternion.AngleAxis instead. This function was deprecated because it uses radians instead of degrees")]
	public void SetAxisAngle(Vector3d axis, double angle)
	{
		this = Quaternion.AxisAngle(axis, angle);
	}
	*/
		#endregion
		public override int GetHashCode()
		{
			return x.GetHashCode() ^ y.GetHashCode() << 2 ^ z.GetHashCode() >> 2 ^ w.GetHashCode() >> 1;
		}
		public override bool Equals(object other)
		{
			if (!(other is QuaternionD))
			{
				return false;
			}
			var quaternion = (QuaternionD)other;
			return x.Equals(quaternion.x) && y.Equals(quaternion.y) && z.Equals(quaternion.z) && w.Equals(quaternion.w);
		}
		public bool Equals(QuaternionD other)
		{
			return x.Equals(other.x) && y.Equals(other.y) && z.Equals(other.z) && w.Equals(other.w);
		}
		public static QuaternionD operator *(QuaternionD lhs, QuaternionD rhs)
		{
			return new QuaternionD(lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y, lhs.w * rhs.y + lhs.y * rhs.w + lhs.z * rhs.x - lhs.x * rhs.z, lhs.w * rhs.z + lhs.z * rhs.w + lhs.x * rhs.y - lhs.y * rhs.x, lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z);
		}
		public static Vector3d operator *(QuaternionD rotation, Vector3d point)
		{
			var num = rotation.x * 2;
			var num2 = rotation.y * 2;
			var num3 = rotation.z * 2;
			var num4 = rotation.x * num;
			var num5 = rotation.y * num2;
			var num6 = rotation.z * num3;
			var num7 = rotation.x * num2;
			var num8 = rotation.x * num3;
			var num9 = rotation.y * num3;
			var num10 = rotation.w * num;
			var num11 = rotation.w * num2;
			var num12 = rotation.w * num3;
			Vector3d result;
			result.x = (1 - (num5 + num6)) * point.x + (num7 - num12) * point.y + (num8 + num11) * point.z;
			result.y = (num7 + num12) * point.x + (1 - (num4 + num6)) * point.y + (num9 - num10) * point.z;
			result.z = (num8 - num11) * point.x + (num9 + num10) * point.y + (1 - (num4 + num5)) * point.z;
			return result;
		}
		public static bool operator ==(QuaternionD lhs, QuaternionD rhs)
		{
			return Dot(lhs, rhs) > 0.999999;
		}
		public static bool operator !=(QuaternionD lhs, QuaternionD rhs)
		{
			return Dot(lhs, rhs) <= 0.999999;
		}
		#region Implicit conversions to and from Unity's Quaternion
		public static implicit operator Quaternion(QuaternionD me)
		{
			return new Quaternion((float)me.x, (float)me.y, (float)me.z, (float)me.w);
		}
		public static implicit operator QuaternionD(Quaternion other)
		{
			return new QuaternionD(other.x, other.y, other.z, other.w);
		}
		#endregion
	}
}