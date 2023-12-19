using System;
using Avionic.Math;
using UnityEngine;

namespace Avionic.Math
{
    [Serializable]
    public class PositionD : MonoBehaviour
    {
        public Vector3d EditorPosition;
        private Transform _target;
        
        public UniverseContainer Universe;
        public Vector3d Value
        {
            get => _value;
            set
            {
                _value = value;
                UpdatePosition();
            }
        }

        private Vector3d _value;

        private void Awake()
        {
            _target = transform;
            
            if (Universe == null)
            {
                Universe = GetComponentInParent<UniverseContainer>();
                
                if (Universe == null)
                {
                    enabled = false;
                    Debug.LogError("[PositionD] This component requires to be within a universe container");
                    return;
                }
            }

            Value = EditorPosition;

            Universe.OffsetChanged += _ => UpdatePosition();
            Universe.ScaleChanged  += _ => UpdatePosition();
            UpdatePosition();
        }

        public Quaternion Rotation => _target.rotation;

        public void MatchRotation(PositionD other)
        {
            _target.rotation = other.Rotation;
        }

        private void Update()
        {
            // Editor update:
            EditorPosition = transform.position - Universe.Offset;
        }

        private void UpdatePosition()
        {
            if (_target == null) return;
            _target.position = (Value + Universe.Offset) * Universe.Scale;
        }
    }
}