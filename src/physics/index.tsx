import { ReactNode } from 'react'
import { useFrame } from '@react-three/fiber'
import { useSimulationStore } from '../state/simulationStore'

// TODO: Refactor to use react-three-rapier Physics component API
// World and RigidBodyDesc are not exported from @react-three/rapier v2.x
// Physics should be managed through the Physics component

export const usePhysicsBody = (descriptor: any) => {
  // TODO: Implement using react-three-rapier API
  // Using any for now since RigidBody type isn't easily accessible
  return { ref: undefined as unknown as React.RefObject<any> }
}

export const PhysicsProvider = ({ children, mode }: { children: ReactNode; mode: string }) => {
  const isPlaying = useSimulationStore((state) => state.isPlaying)
  const simSpeed = useSimulationStore((state) => state.simSpeed)
  const time = useSimulationStore((state) => state.time)
  const setTime = useSimulationStore((state) => state.setTime)

  useFrame((_, delta) => {
    if (isPlaying) {
      // TODO: Use react-three-rapier Physics API for stepping
      setTime(time + delta * simSpeed)
    }
  })

  return <>{children}</>
}

export const createGliderBody = async () => {
  // TODO: Implement using react-three-rapier API
  return null
}

export const createSectionBody = async (position: { x: number; y: number; z: number }) => {
  // TODO: Implement using react-three-rapier API
  return null
}

export const stepPhysics = async (deltaTime: number) => {
  // TODO: Implement using react-three-rapier API
}
