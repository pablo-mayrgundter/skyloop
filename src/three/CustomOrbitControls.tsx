import { useThree, useFrame } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import { useRef, useMemo, useEffect } from 'react'
import { OrbitControls as OrbitControlsImpl } from 'three-stdlib'
import { Vector3, Box3, Sphere } from 'three'
import { EARTH_RADIUS_KM, LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle, LOOP_RADIUS_KM } from '../config/simulation'

// Calculate zoom speed based on distance from camera to loop target (in meters)
// Exponential slowdown as approaching loop, exponential speedup as zooming out
const calculateZoomSpeed = (distance: number): number => {
  // Distance from camera to loop target (in km)
  // Small distance = close to loop (slow zoom), large distance = far from loop (fast zoom)
  const clampedDistance = Math.max(0.0001, distance) // Minimum 0.0001 km = 0.1 meters
  
  // Logarithmic scaling: zoom speed decreases exponentially as we get closer
  // At far distances (100km): speed = 2.0 (fast zoom out)
  // At close distances (0.001km = 1m): speed ≈ 0.0001 (very slow zoom in, ~0.1 m/s)
  const normalizedDistance = Math.log10(clampedDistance)
  const minSpeed = 0.0001  // Very slow when close to loop (~0.1 m/s)
  const maxSpeed = 2.0     // Fast when far from loop
  const speedRange = maxSpeed - minSpeed
  
  // Map log distance to speed (closer = slower, farther = faster)
  // log10(100) ≈ 2, log10(0.001) ≈ -3
  // Scale from -3 to 2, map to minSpeed to maxSpeed
  const speed = minSpeed + speedRange * ((normalizedDistance + 3) / 5)
  
  return Math.max(minSpeed, Math.min(maxSpeed, speed))
}

export const CustomOrbitControls = (props: React.ComponentProps<typeof OrbitControls>) => {
  const controlsRef = useRef<OrbitControlsImpl>(null)
  const { camera, scene } = useThree()
  const hasFittedRef = useRef(false)
  
  // Target a point on Earth's surface near the loop (for viewing components)
  // This gives a better view of the loop components on the surface
  const surfaceTarget = useMemo(() => {
    return interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, EARTH_RADIUS_KM)
  }, [])
  
  // Target Earth center (0, 0, 0) as fallback
  const earthTarget = useMemo(() => new Vector3(0, 0, 0), [])
  
  // Target a point on the loop (midpoint of the great circle path)
  const loopTarget = useMemo(() => {
    return interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, LOOP_RADIUS_KM)
  }, [])
  
  // Fit camera to entire Earth on initialization
  useEffect(() => {
    // Use a small delay to ensure controls are fully initialized
    const timer = setTimeout(() => {
      const controls = controlsRef.current
      if (!controls || hasFittedRef.current) return
      
      // OrbitControlsImpl extends CameraControls, so it should have fitToSphere/fitToBox
      const cameraControls = controls as any
      
      // Create bounding sphere for Earth (centered at origin)
      const earthSphere = new Sphere(new Vector3(0, 0, 0), EARTH_RADIUS_KM)
      
      // Try fitToSphere first (more appropriate for spherical Earth)
      if (cameraControls && typeof cameraControls.fitToSphere === 'function') {
        cameraControls.fitToSphere(earthSphere, false)
        // After fitting, adjust target to look at surface components
        controls.target.copy(surfaceTarget)
        controls.update()
        hasFittedRef.current = true
        return
      }
      
      // Fallback to fitToBox
      const earthBox = new Box3(
        new Vector3(-EARTH_RADIUS_KM, -EARTH_RADIUS_KM, -EARTH_RADIUS_KM),
        new Vector3(EARTH_RADIUS_KM, EARTH_RADIUS_KM, EARTH_RADIUS_KM)
      )
      
      if (cameraControls && typeof cameraControls.fitToBox === 'function') {
        cameraControls.fitToBox(earthBox, false, {
          paddingLeft: 0.1,
          paddingRight: 0.1,
          paddingTop: 0.1,
          paddingBottom: 0.1
        })
        // After fitting, adjust target to look at surface components
        controls.target.copy(surfaceTarget)
        controls.update()
        hasFittedRef.current = true
        return
      }
      
      // Final fallback: manually calculate camera position to frame Earth
      // Position camera at 2-3x Earth radius (in space, looking at Earth)
      const cameraDistance = EARTH_RADIUS_KM * 2.5 // 2.5x Earth radius (~15,927km)
      
      // Position camera at good viewing angle (slightly elevated, angled view)
      camera.position.set(
        cameraDistance * 0.7,  // X: offset for angled view
        cameraDistance * 0.5,  // Y: elevated view  
        cameraDistance * 0.7   // Z: offset for angled view
      )
      
      // Look at Earth surface near the loop (better view of components)
      camera.lookAt(surfaceTarget)
      controls.target.copy(surfaceTarget)
      controls.update()
      hasFittedRef.current = true
    }, 100) // Small delay to ensure controls are ready
    
    return () => clearTimeout(timer)
  }, [camera, earthTarget])
  
  useFrame(() => {
    const controls = controlsRef.current
    if (!controls) return
    
    // After initial fit, allow targeting loop elements
    // But keep Earth center as default target for better overall view
    if (!hasFittedRef.current) return
    
    // Prevent camera from going through Earth surface
    const cameraDistanceFromOrigin = camera.position.length()
    if (cameraDistanceFromOrigin < EARTH_RADIUS_KM + 0.01) { // 0.01 km = 10 meters above surface
      // Push camera back to just above surface
      const direction = camera.position.clone().normalize()
      camera.position.copy(direction.multiplyScalar(EARTH_RADIUS_KM + 0.01))
      controls.update()
    }
    
    // Calculate distance from camera to target (Earth center or loop)
    const currentTarget = controls.target.distanceTo(loopTarget) < 1 ? loopTarget : earthTarget // 1 km tolerance
    const distance = camera.position.distanceTo(currentTarget)
    
    // Update zoom speed dynamically based on distance
    // Speed slows down exponentially as approaching, speeds up as zooming out
    const zoomSpeed = calculateZoomSpeed(distance)
    controls.zoomSpeed = zoomSpeed
  })
  
  return <OrbitControls ref={controlsRef} target={surfaceTarget.toArray()} {...props} />
}

