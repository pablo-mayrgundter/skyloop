export type SimulationMode = 'firstPerson' | 'birdseye'

export type SkyloopModuleId =
  | 'main'
  | 'sections'
  | 'tethers'
  | 'tube'
  | 'gliders'
  | 'atmosphere'
  | 'weather'
  | 'dayNight'
  | 'bom'
  | 'economics'

export type OverlayFlags = {
  showAtmosphere: boolean
  showWind: boolean
  showTethers: boolean
  showSections: boolean
  showGliders: boolean
  showDayNight: boolean
}

export type AtmosphericLayer = {
  name: string
  altitudeKm: number
  density: number
}

export type WindProfile = {
  altitudeKm: number
  speed: number
  directionDeg: number
}

export type CostBreakdown = {
  category: string
  costMillions: number
}
