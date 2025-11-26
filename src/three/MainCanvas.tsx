import { Canvas } from '@react-three/fiber'
import { useMemo } from 'react'
import SceneRoot from './SceneRoot'
import { PhysicsProvider } from '../physics'
import { useSimulationStore } from '../state/simulationStore'
import { CustomOrbitControls } from './CustomOrbitControls'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO, LOOP_RADIUS_KM, EARTH_RADIUS_KM } from '../config/simulation'
import { Vector3 } from 'three'

// Zoom distances in km
// Min distance: 0.001 km = 1 meter from loop element (human scale viewing)
// Max distance: 3x Earth radius for far space view
const HUMAN_SCALE_DISTANCE = 0.001 // 0.001 km = 1 meter from loop element
const MAX_ZOOM_DISTANCE = EARTH_RADIUS_KM * 3 // 3x Earth radius (~19,113km)

// Camera near and far planes (in km)
// Near: 0.0001 km = 0.1m for close-up viewing (human scale)
// Far: 20,000km to see entire Earth and beyond
const CAMERA_NEAR = 0.0001
const CAMERA_FAR = 20000 // 20,000km

// Calculate initial camera position for birds-eye view
// Position camera 2-3x Earth radius away in space, looking at Earth surface/components
const getInitialCameraPosition = (): [number, number, number] => {
  // Position camera at 2.5x Earth radius (about 15,927km from center)
  const cameraDistance = EARTH_RADIUS_KM * 2.5
  
  // Position camera looking at the loop area on Earth's surface
  // Use a good viewing angle (slightly above equator, looking down)
  // This gives a nice perspective view of Earth and the loop components
  const cameraPos = new Vector3(
    cameraDistance * 0.7,  // X: offset for angled view
    cameraDistance * 0.5,  // Y: elevated view
    cameraDistance * 0.7   // Z: offset for angled view
  )
  
  return [cameraPos.x, cameraPos.y, cameraPos.z]
}

const MainCanvas = () => {
  const mode = useSimulationStore((state) => state.mode)
  
  // Calculate initial camera position based on loop location
  const initialCameraPos = useMemo(() => getInitialCameraPosition(), [])

  return (
    <Canvas 
      shadows 
      camera={{ 
        position: initialCameraPos,
        fov: 55,
        near: CAMERA_NEAR,
        far: CAMERA_FAR
      }}
    >
      <color attach="background" args={[0.02, 0.04, 0.09]} />
      <ambientLight intensity={0.2} />
      {/* Position light relative to Earth scale */}
      <directionalLight position={[EARTH_RADIUS_KM * 2, EARTH_RADIUS_KM * 2, EARTH_RADIUS_KM]} intensity={1.2} />
      <PhysicsProvider mode={mode}>
        <SceneRoot />
      </PhysicsProvider>
      {/* Only enable OrbitControls in birds-eye mode - other modes have their own camera control */}
      {mode === 'birdseye' && (
        <CustomOrbitControls 
          enablePan={false} 
          minDistance={HUMAN_SCALE_DISTANCE} 
          maxDistance={MAX_ZOOM_DISTANCE}
          // Prevent going through Earth surface
          minPolarAngle={0}
          maxPolarAngle={Math.PI}
        />
      )}
      {/* StatsGl not available in drei v10 - can add stats.js separately if needed */}
    </Canvas>
  )
}

export default MainCanvas
