import { useFrame, useThree } from '@react-three/fiber'
import { useRef } from 'react'
import { PerspectiveCamera, Vector3 } from 'three'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO, LOOP_RADIUS_KM } from '../config/simulation'
import { useSimulationStore } from '../state/simulationStore'

export const useFirstPersonCamera = () => {
  const camera = useThree((state) => state.camera as PerspectiveCamera)
  const time = useSimulationStore((state) => state.time)
  const simSpeed = useSimulationStore((state) => state.simSpeed)
  const lookAheadRef = useRef(new Vector3())
  
  useFrame(() => {
    const loopTime = (time * simSpeed) % 300
    const t = (loopTime % 60) / 60
    
    // Current position along the loop
    const currentPos = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t, LOOP_RADIUS_KM)
    
    // Look ahead along the path (slightly ahead for smooth camera movement)
    const lookAheadT = Math.min(1, t + 0.01)
    const lookAheadPos = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, lookAheadT, LOOP_RADIUS_KM)
    
    // Set camera position
    camera.position.copy(currentPos)
    
    // Look ahead along the path
    camera.lookAt(lookAheadPos)
    
    // Store look ahead for potential use
    lookAheadRef.current.copy(lookAheadPos)
  })
}

export const FirstPersonRig = () => {
  useFirstPersonCamera()
  return null
}

export const useFollowCamera = () => {
  const camera = useThree((state) => state.camera as PerspectiveCamera)
  const time = useSimulationStore((state) => state.time)
  const simSpeed = useSimulationStore((state) => state.simSpeed)
  const pathRef = useRef<Vector3[] | null>(null)
  
  // Pre-calculate path once
  if (!pathRef.current) {
    const path: Vector3[] = []
    const segments = 64
    for (let i = 0; i <= segments; i += 1) {
      path.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, i / segments, LOOP_RADIUS_KM))
    }
    pathRef.current = path
  }
  
  useFrame(() => {
    const path = pathRef.current!
    
    // Calculate glider position (same logic as LoopGliders)
    const loopTime = (time * simSpeed) % path.length
    const idx = Math.floor(loopTime)
    const t = loopTime - idx
    const current = path[idx % path.length]
    const next = path[(idx + 1) % path.length]
    const gliderPosition = current.clone().lerp(next, t)
    
    // Calculate glider forward direction (normalized)
    const gliderForward = next.clone().sub(current).normalize()
    
    // Calculate camera offset: behind and slightly above the glider
    // Distance: 0.1 km (100m) behind, 0.05 km (50m) above
    const followDistance = 0.1 // km behind glider
    const followHeight = 0.05 // km above glider
    
    // Calculate right vector (perpendicular to forward and up)
    const worldUp = new Vector3(0, 1, 0)
    const right = new Vector3().crossVectors(gliderForward, worldUp).normalize()
    
    // Calculate actual up vector (perpendicular to forward and right)
    const up = new Vector3().crossVectors(right, gliderForward).normalize()
    
    // Position camera: glider position + offset behind and up
    const cameraOffset = gliderForward
      .clone()
      .multiplyScalar(-followDistance) // Behind
      .add(up.clone().multiplyScalar(followHeight)) // Above
    
    camera.position.copy(gliderPosition).add(cameraOffset)
    
    // Look at glider position (or slightly ahead for better view)
    const lookAhead = gliderForward.multiplyScalar(0.05) // 50m ahead
    const lookTarget = gliderPosition.clone().add(lookAhead)
    camera.lookAt(lookTarget)
  })
}

export const FollowRig = () => {
  useFollowCamera()
  return null
}

export const BirdsEyeRig = ({ overlays }: { overlays: Record<string, boolean> }) => {
  // Don't manipulate camera - let OrbitControls handle it
  // This component just renders overlay elements
  return (
    <group>
      {overlays.showWind && (
        <mesh position={[0, 1.6, 0]}>
          <coneGeometry args={[0.05, 0.2, 6]} />
          <meshStandardMaterial color="#74f0ed" transparent opacity={0.6} />
        </mesh>
      )}
    </group>
  )
}
